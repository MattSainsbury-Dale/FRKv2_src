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
## knots_squared     : spBayes
if (quick) {
  ARGS <- list(max.edge.interior = 50, nres = 1, k = 20, n.neighbours = 2, knots_squared = 9) 
} else {
  ARGS <- list(max.edge.interior = 4, nres = 4, k = 2250, n.neighbours = 15, knots_squared = 20^2)
}

quiet_source <- function(...) suppressMessages(source(...))

## Load the model fitting and prediction functions
mapply(quiet_source, paste0("./scripts/MODIS_modelling_fns/", PACKAGES, ".R"))

# suppressMessages(
#   mapply(source, paste0("./scripts/MODIS_modelling_fns/", PACKAGES, ".R"))
# )

## Load user defined functions 
## These include functions to compute diagnostics, convert data frames from wide
## to long format in a convenient way for this study, and some other 
## miscellaneous functions to make the comparison study script easier to read.
source("./scripts/MODIS_diagnostic_and_misc_fns.R")
source("./scripts/Plotting_helpers/Plotting_helpers.R")
source("./scripts/Plotting_helpers/MODIS_plotting_fns.R")


# ---- Load MODIS data  ----

data("MODIS_cloud_df") 
df <- MODIS_cloud_df 
rm(MODIS_cloud_df)


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
missing_fwk <- "MAR"
for (i in 1:number_of_runs_MAR) {
  seed <- i # Controls the training set
  source("./scripts/MODIS_analysis.R")
  
  df_train$Run <- i
  df_test$Run <- i
  df_train$Sampling_scheme <- "MAR"
  df_test$Sampling_scheme <- "MAR"
  
  MAR_df_train_list[[i]] <- df_train
  MAR_df_test_list[[i]] <- df_test
  MAR_times[i, ] <- c(times[PACKAGES], "MAR", i)
}

## "Missing at block"
block_list <- block_df_train_list <- block_df_test_list <- list()
block_times <- data.frame(matrix(ncol = length(PACKAGES) + 2, nrow = 0))
colnames(block_times) <- c(PACKAGES, "Sampling_scheme", "Run")
missing_fwk <- "block"
for (i in 1:number_of_runs_block) {
  seed <- i + 1 # NB: i + 1 gives an interesting block position when number_of_runs_block == 1
  source("./scripts/MODIS_analysis.R")
  
  df_train$Run <- i
  df_test$Run <- i
  df_train$Sampling_scheme <- "block"
  df_test$Sampling_scheme <- "block"
  
  block_df_train_list[[i]] <- df_train
  block_df_test_list[[i]] <- df_test
  block_times[i, ] <- c(times[PACKAGES], "block", i)
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

