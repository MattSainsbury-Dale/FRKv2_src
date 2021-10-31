source("./scripts/Utility_fns.R")


## Use very-low-dimensional representations of the models to establish that the code works? 
quick <- check_quick()

## Packages used (use whichever subset you please)
PACKAGES <- c(
  "FRK",
  "INLA",
  "mgcv",
  "spNNGP",
  "spBayes"
)

# ---- Load packages and user-defined functions ----


suppressMessages({
## Packages used for scoring rules and plotting
library("pROC")
library("ggplot2")
library("ggpubr")
library("dplyr")
library("sp")
library("stringr")
library("tidyr")
  
## Packages required by INLA
library("foreach")
library("splancs")
library("rgdal")

## Packages used in the comparison study
library("FRK") 
library("mgcv")
library("INLA")
library("spNNGP")
library("spBayes")
})


options(dplyr.summarise.inform = FALSE) # Suppress summarise info

## Arguments to packages: These arguments in some way define the number of 
## basis functions used for each model. For instance, nres specifies the number of 
## basis-functions resolutions to use with FRK(), and knot_squared specifies the
## number of knots to use for spBayes. 
## max.edge.interior : INLA
## nres              : FRK
## k                 : mgcv
## n.neighbours      : spNNGP
## knots             : spBayes (should be a square number)
if (quick) {
  ARGS <- list(max.edge.interior = 50, nres = 1, k = 20, n.neighbours = 2, knots = 3^2) 
} else {
  ARGS <- list(max.edge.interior = 4, nres = 4, k = 2250, n.neighbours = 15, knots = 20^2)
}



## Load the model fitting and prediction functions
## (assign to dummy variable to prevent output)
dummy <- mapply(source, paste0("./scripts/MODIS_modelling_fns/", PACKAGES, ".R"))



# ---- Helper functions ----

## These include functions to compute diagnostics, convert data frames from wide
## to long format in a convenient way for this study, and some plotting 
## helper functions.

## Define the functions we will use for comparing across methods
BrierScore <- function(Z, pred) mean((Z - pred)^2)
AUC <- function(Z, pred) suppressMessages(as.numeric(pROC::auc(Z, pred)))

compute_diagnostics_MODIS <- function(df) {
  summarise(df, Brier = BrierScore(z, pred), AUC = AUC(z, pred))
}

## Creates long form dataframe, useful for diagnostics and plotting
long_prediction_data_frame <- function(df) {
  data.frame(
    Method = rep(PACKAGES, each = nrow(df)),
    x    = rep(df$x, times = length(PACKAGES)),
    y    = rep(df$y, times = length(PACKAGES)),
    pred   = c(as.matrix(df[, paste0("pred_", PACKAGES)]))
  )
}

## Define a function needed for the missing-block sampling scheme.
## Samples a rectangular Block of width w and height h from a rectangular spatial 
## region defined by the coordinates x and y in df. 
## It randomly samples a point within D which acts as a 
## lower-left vertex for the Block. The sample region for this vertex depends on 
## the size of the Block; we cannot have it too close to the right- or top-edge of 
## D, otherwise the Block will extend beyond the considered boundaries of D. 
sampleBlock <- function(df, w, h) {
  ## Retain only those coordinates which would result in the entire Block being inside D
  df <- df %>% subset(x < (max(x) - w) & y < (max(y) - h))
  
  ## Sample a point from the valid domain which will act as the lower-left vertex
  v1 <- df[sample(1:nrow(df), 1), c("x", "y")]
  
  ## Now create the upper-right vertex
  v2 <- v1 + c(w, h)
  
  return(c(xmin = v1[, 1], 
           ymin = v1[, 2], 
           xmax = v2[, 1], 
           ymax = v2[, 2]))
}


cloud_colour    = "orange"  # "white"
no_cloud_colour = "blue"    # "black"
missing_colour  = "white"   # "#BFD5E3"
midpoint_colour = "gray" 

change_font_size <- function(gg) {
  gg + theme(axis.text = element_text(size = 10),
             axis.title = element_text(size = 13), 
             legend.text = element_text(size = 10), 
             strip.text = element_text(size = 12))
}

training_data_background <- theme(
  panel.background = element_rect(fill = "white", colour = "white"), # "#6D9EC1"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

discrete_cloud_theme <- theme(legend.key=element_rect(
  # fill = "gray", colour = "gray", # use this if you colour the clouds white
  fill = "white", colour = "white", 
  size = 2))


discrete_cloud_scale <- scale_fill_gradient(low = no_cloud_colour, 
                                            high = cloud_colour,
                                            breaks = c(0, 1),
                                            guide = "legend",
                                            labels = c("No Cloud", "Cloud"),
                                            name = "")


common_layers <- ggplot() + theme_bw() + coord_fixed() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + 
  labs(x = bquote(s[1]), y = bquote(s[2]))



# ---- Load MODIS data  ----

# FIXME: Just keep MODIS_cloud_df name
data("MODIS_cloud_df") 
df <- MODIS_cloud_df 
rm(MODIS_cloud_df)


# ---- Analysis function ----

MODIS_analysis <- function(df,  Sampling_scheme, PACKAGES, ARGS) {
  
  # ---- Create training and test sets: df_train and df_test ----
  
  ## First compute observed indices.
  ## The way we sample the observed indices depends on whether the sampling scheme
  ## is missing at random, or "missing at block".
  if(!exists("seed")) seed <- 1 
  RNGversion("3.6.0")
  set.seed(seed)
  
  if(!exists("Sampling_scheme")) {
    Sampling_scheme <- "MAR"
  }
  if(Sampling_scheme == "MAR") {
    n <- 6000
    obs.id <- sample(1:nrow(df), n, replace = FALSE) 
  } else if (Sampling_scheme == "block") {
    block  <- sampleBlock(df, w = 30, h = 30)
    obs.id <- which(!with(df, (x > block["xmin"] & x <= block["xmax"]) & (y > block["ymin"] & y <= block["ymax"])))
  }
  
  ## Based on the observed indices, construct the unobserved indices, 
  ## training set, and test set.
  unobs.id <- (1:nrow(df))[-obs.id] 
  df_train <- df[obs.id, ]   
  df_test  <- df[unobs.id, ]
  
  # ---- Run the models ----
  
  times         <- list()         ## store fitting times
  model_objects <- list()  ## store model objects (which may be of use later)
  
  ## More informative information about the sampling framework
  if(Sampling_scheme == "MAR") {
    msg <- "missing at random (MAR)"
  } else if (Sampling_scheme == "block") {
    msg <- "missing block (MB)"
  }
  
  cat(paste("\n#### Starting the MODIS comparison using the"), msg, "sampling scheme ####\n")
  
  print_start_msg <- function(pkg) cat(paste("\n---- Running", pkg, "on", msg, "MODIS data ----\n\n"))
  
  print_run_time_MODIS <- function(pkg, times) {
    cat(paste("Finished", pkg, "analysis in", round(times[[pkg]]["elapsed"] / 60, 4), "minutes.\n"))
  }
  
  if("FRK" %in% PACKAGES) {
    pkg <- "FRK"
    print_start_msg(pkg)
    times$FRK <- system.time({
      model_objects$FRK <- MODIS_FRK_fit(df_train, nres = ARGS$nres)
      df_test$pred_FRK  <- MODIS_FRK_pred(df_test, model_objects$FRK)
    })
    print_run_time_MODIS(pkg, times)
  }
  
  
  if("INLA" %in% PACKAGES) {
    pkg <- "INLA"
    print_start_msg(pkg)
    times$INLA <- system.time(
      df_test$pred_INLA <- MODIS_INLA(pred_locs = df_test, df_train = df_train, 
                                      max.edge.interior = ARGS$max.edge.interior)
    )
    print_run_time_MODIS(pkg, times)
  }
  
  if("mgcv" %in% PACKAGES) {
    pkg <- "mgcv"
    print_start_msg(pkg)
    times$mgcv <- system.time({
      model_objects$mgcv <- MODIS_mgcv_fit(df_train, k = ARGS$k)
      df_test$pred_mgcv  <- MODIS_mgcv_pred(df_test, model_objects$mgcv)
    })
    print_run_time_MODIS(pkg, times)
  }
  
  
  if("spNNGP" %in% PACKAGES) {
    pkg <- "spNNGP"
    print_start_msg(pkg)
    times$spNNGP <- system.time({
      model_objects$spNNGP <- MODIS_spNNGP_fit(df_train, n.neighbours = ARGS$n.neighbours)
      df_test$pred_spNNGP  <- MODIS_spNNGP_pred(df_test, model_objects$spNNGP)
    })
    print_run_time_MODIS(pkg, times)
  }
  
  
  if("spBayes" %in% PACKAGES) {
    pkg <- "spBayes"
    print_start_msg(pkg)
    times$spBayes <- system.time({
      model_objects$spBayes <- MODIS_spBayes_fit(df_train, knots = ARGS$knots)
      df_test$pred_spBayes <- MODIS_spBayes_pred(df_test, model_objects$spBayes)
    })
    print_run_time_MODIS(pkg, times)
  }
  
  
  times <- sapply(times, function(x) unname(x["elapsed"]))
  times <- times[PACKAGES] 
  
  
  # ---- Fitted values ----
  
  ## Compute fitted values so that we can predict over D
  
  if("INLA" %in% PACKAGES)
    df_train$pred_INLA <- INLA_fitted_values # NB: INLA_fitted_values computed within MODIS_INLA()
  
  if("FRK" %in% PACKAGES)
    df_train$pred_FRK <- FRK_fitted_values # NB: FRK_fitted_values computed within MODIS_FRK_pred()
  
  if("mgcv" %in% PACKAGES)
    df_train$pred_mgcv <- model_objects$mgcv$fitted.values
  
  if("spNNGP" %in% PACKAGES)
    df_train$pred_spNNGP <- rowMeans(plogis(model_objects$spNNGP$y.hat.samples))
  
  if("spBayes" %in% PACKAGES)
    df_train$pred_spBayes <- fitted_values_spBayes(model_objects$spBayes)
  
  
  # ---- Output ----
  
  ## The controlling script now uses df_train and df_test, which now also contains
  ## fitted values and predictions, to compute diagnostics and ROC curves (using 
  ## df_test) and prediction maps over D (using rbind(df_train, df_test))
  
  results <- list(df_test = df_test, df_train = df_train, times = times)
  
  if (Sampling_scheme == "block") results$blocks <- blocks
  
  return(results)
}

# ---- Run comparison study and save results ----

## Run the study with multiple training data sets. We source MODIS_anlaysis.R in a loop, 
## changing the seed with which the training data is defined in each iteration. 
## A list of diagnostic scores and plots are recorded at each iteration.
## In the paper, we only did one run for each sampling scheme because 
## i) the results didn't change very much between runs, 
## ii) the run-time is already very large even with only 1 run, and
## iii) the exposition is complicated by multiple runs.
## We keep the code here for future reference, and if anyone wants to test the 
## results using multiple runs. 

number_of_runs_MAR <- 1
number_of_runs_block <- 1



## Missing at random
MAR_df_train_list <- MAR_df_test_list <- list()
MAR_times <- data.frame(matrix(ncol = length(PACKAGES) + 2, nrow = 0))
colnames(MAR_times) <- c(PACKAGES, "Sampling_scheme", "Run")
Sampling_scheme <- "MAR"
for (i in 1:number_of_runs_MAR) {
  seed <- i # Controls the training set
  source("./scripts/MODIS_analysis.R")
  
  df_train$Run <- i
  df_test$Run <- i
  df_train$Sampling_scheme <- Sampling_scheme
  df_test$Sampling_scheme <- Sampling_scheme
  
  MAR_df_train_list[[i]] <- df_train
  MAR_df_test_list[[i]] <- df_test
  MAR_times[i, ] <- c(times[PACKAGES], Sampling_scheme, i)
}

## "Missing at block"
block_list <- block_df_train_list <- block_df_test_list <- list()
block_times <- data.frame(matrix(ncol = length(PACKAGES) + 2, nrow = 0))
colnames(block_times) <- c(PACKAGES, "Sampling_scheme", "Run")
Sampling_scheme <- "block"
for (i in 1:number_of_runs_block) {
  seed <- i + 1 # NB: i + 1 gives an interesting block position when number_of_runs_block == 1
  source("./scripts/MODIS_analysis.R")
  
  df_train$Run <- i
  df_test$Run <- i
  df_train$Sampling_scheme <- Sampling_scheme
  df_test$Sampling_scheme <- Sampling_scheme
  
  block_df_train_list[[i]] <- df_train
  block_df_test_list[[i]] <- df_test
  block_times[i, ] <- c(times[PACKAGES], Sampling_scheme, i)
  block_list[[i]] <- block # record the block for plotting purposes
}

blocks <- do.call("rbind", block_list) %>% 
  as.data.frame()  %>%
  mutate(Run = 1:length(block_list))
write.csv(blocks,
          file = "./intermediates/MODIS_blocks.csv",
          row.names = FALSE)

## Combine data frames, and save as a .csv
tmp1 <- do.call("rbind", MAR_df_train_list)
tmp2 <- do.call("rbind", block_df_train_list)
all_df_train <- rbind(tmp1, tmp2)
write.csv(all_df_train,
          file = "./intermediates/MODIS_all_df_train.csv",
          row.names = FALSE)

tmp1 <- do.call("rbind", MAR_df_test_list)
tmp2 <- do.call("rbind", block_df_test_list)
all_df_test <- rbind(tmp1, tmp2)
write.csv(all_df_test,
          file = "./intermediates/MODIS_all_df_test.csv",
          row.names = FALSE)

times <- rbind(MAR_times, block_times)
times <- times %>% gather(Method, time, all_of(PACKAGES))
times$Run <- as.integer(times$Run)
times$time <- as.numeric(times$time)
times$time <- times$time/60 # convert to minutes
write.csv(times,
          file = "./intermediates/times.csv",
          row.names = FALSE)

## ---- Re-load results (useful to change plots without running models) ----

blocks <- read.csv("./intermediates/MODIS_blocks.csv")
all_df_train <- read.csv("./intermediates/MODIS_all_df_train.csv")
all_df_test <- read.csv(file = "./intermediates/MODIS_all_df_test.csv")
times <- read.csv(file = "./intermediates/times.csv")


## ---- Training and testing data visualisation ----

## See MODIS_plotting_fns.R for the cloud colour definitions

## Plot the unthresholded version of the data
g_original_data <- common_layers +
  geom_raster(data = df, aes(x, y, fill = z_unthresholded)) +
  scale_fill_gradient(low = no_cloud_colour, high = cloud_colour) +
  labs(fill = "Radiance")

## Plot the thresholded data (which we use for our analyses)

## Full data set
g_thresholded_data <- common_layers +
  geom_raster(data = df, aes(x, y, fill = z)) +
  discrete_cloud_scale + discrete_cloud_theme

## Plot the training data
## First create a data frame which facilitates faceting
tmp1 <- all_df_train %>% mutate(set = "training")
tmp2 <- all_df_test %>% mutate(set = "test", time = NULL)
tmp <- rbind(tmp1, tmp2) %>%
  mutate(facet_var = paste(Sampling_scheme, set)) %>%
  mutate(facet_var = factor(
    facet_var, 
    levels = c("MAR training", "MAR test", "block training", "block test"),
    labels = c("MR: training set", "MR: test set", "MB: training set", "MB: test set")
    )) %>%
  filter(Run == 1)

figure <- (common_layers +
  geom_raster(data = tmp, aes(x = x, y = y, fill = z))  +
  facet_wrap(~facet_var, nrow = 1)) %>% 
  change_font_size + 
  discrete_cloud_scale + discrete_cloud_theme + training_data_background

ggsave(figure,
       path = "./results/", filename = "4_1_MODIS_data.png", device = "png", 
       width = 10, height = 2.2)

# ---- Compute diagnostics (Brier score and AUC) ----

## Convert all_df_test to long format, and add run times
all_df_test <- all_df_test %>% 
  gather(Method, pred, paste0("pred_", PACKAGES)) %>%
  mutate(Method = str_remove(Method, "pred_")) %>% 
  left_join(times, by = c("Run", "Sampling_scheme", "Method"))

## Compute diagnostics, split by run, method, and sampling scheme
diagnostics <- all_df_test %>%
  dplyr::group_by(Sampling_scheme, Run, Method) %>%
  dplyr::summarise(Brier = BrierScore(z, pred),
            AUC = AUC(z, pred), 
            time = time[1])

## Compute average diagnostics, split by method and sampling scheme, but averaged over all runs
## (the averaging is redundant for the paper, because we do only one run)
diagnostics <- diagnostics %>%  
  dplyr::group_by(Sampling_scheme, Method) %>%
  dplyr::summarise(
    Brier = mean(Brier),
    AUC = mean(AUC), 
    time = mean(time)) %>% 
  as.data.frame()

write.csv(diagnostics, row.names = FALSE, 
          file = "./results/4_1_MODIS_diagnostics.csv")

save_html_table(
  diagnostics,
  file = "results/4_1_MODIS_diagnostics.html", 
  caption = "MODIS comparison study"
  )

# ---- Prediction maps ----

## Convert all_df_train to long format
all_df_train <- all_df_train %>% 
  gather(Method, pred, paste0("pred_", PACKAGES)) %>%
  mutate(Method = str_remove(Method, "pred_")) 

## Replace "FRK" with "FRK v2"
all_df_train$Method[all_df_train$Method == "FRK"] <- "FRK v2"
all_df_test$Method[all_df_test$Method == "FRK"] <- "FRK v2"

## Add "Data" as a level of the variable Method, so that the true, 
## withheld data are included in the plots. 
N_test_locations <- nrow(all_df_test) / length(PACKAGES)
tmp <- all_df_test[1:N_test_locations, ]
tmp$Method <- "Data"
tmp$pred <- tmp$z
all_df_test <- rbind(all_df_test, tmp)


plot_predictions <- function(df, scheme) {
  
  if (scheme == "block") {
    facet_type <- facet_grid( ~ Method)
  } else if (scheme == "MAR") {
    facet_type <- facet_wrap( ~ Method, nrow = 2)
  }

  gg <- common_layers %+% filter(df, Sampling_scheme == scheme, Run == 1)
  
  gg <- gg + geom_raster(aes(x, y, fill = pred)) +
    facet_type +
    scale_fill_gradient2(
      midpoint = 0.5, low = no_cloud_colour, 
      mid = midpoint_colour, high = cloud_colour
      ) +
    labs(fill = "")
  
  gg <- gg %>% change_font_size 
  
  return(gg)
}

ggsave(plot_predictions(all_df_test, "MAR"),
       path = "./results/", filename = "4_1_MODIS_MAR_predictions.png", device = "png", 
       width = 10, height = 5.1)

ggsave(plot_predictions(all_df_test, "block"),
       path = "./results/", filename = "4_1_MODIS_block_predictions.png", device = "png",
       width = 12, height = 2.5)   

# ---- ROC curves ----

compute_ROC_objects <- function(df, scheme, run = 1) {
  ROC_list <- list()
  for (method in unique(all_df_test$Method)) {
    ROC_list[[method]] <- all_df_test %>%
      filter(Method == method) %>% 
      filter(Sampling_scheme == scheme, Run == run) %>% 
      pROC::roc(z, pred, quiet = TRUE, auc = FALSE)
  }
  return(ROC_list)
}


block_ROC_list <- compute_ROC_objects(all_df_test, scheme = "block")
MAR_ROC_list   <- compute_ROC_objects(all_df_test, scheme = "MAR")


FRK_idx <- which(PACKAGES == "FRK")

plot_ROC <- function(ROC_list) {
  
  # If present, plot FRK v2 at the end so that it is clearly visible
  if ("FRK v2" %in% names(ROC_list)) ROC_list <- ROC_list[c(PACKAGES[-FRK_idx], "FRK v2")]  
  
  # Remove Data from ROC list (it's just the observed data, so it will have a perfect ROC)
  ROC_list$Data <- NULL
  
  pROC::ggroc(
    ROC_list, 
    legacy.axes = TRUE, size = 0.3) +
    scale_colour_manual(values = 1:5) +
    xlab("False positive rate") + ylab("True positive rate") +
    labs(colour = "Method") +
    geom_abline(slope = 1, intercept = 0, color = "darkgrey", linetype = "dashed") +
    theme_bw() + coord_fixed()
}

block_ROC <- plot_ROC(block_ROC_list)
MAR_ROC   <- plot_ROC(MAR_ROC_list)

figure <- ggpubr::ggarrange(MAR_ROC, block_ROC, 
                    align = "h", common.legend = TRUE, legend = "right")

ggsave(figure, path = "./results/", filename = "4_1_MODIS_ROC.png", device = "png", 
       width = 6, height = 2.6)

