## This script is used to conduct the MODIS comparison study with multiple 
## training data sets. It sources comparison_study.R in a loop, 
## changing the seed with which the training data is defined in each iteration. 
## A list of diagnostic scores and plots are recorded at each iteration.

# ---- Load packages and user-defined functions ----

## Packages used for scoring rules and plotting
library("pROC")
library("ggplot2")
library("ggpubr")
library("dplyr")
library("sp")
library("stringr")
library("tidyr")
library("stargazer")

## Packages used in the comparison study
library("FRK") # FRKTMB branch
library("mgcv")
library("INLA")
library("spNNGP")
library("spBayes")

## Packages used
PACKAGES <- c("FRK", "INLA", "mgcv", "spNNGP", "spBayes")

## Use extremely low-rank versions of the models to establish that the 
## code works (relatively) quickly
LOWRANK <- FALSE
if (LOWRANK) {
  ARGS <- list(max.edge.interior = 50, nres = 1, k = 20, n.neighbours = 2, knots_squared = 9) 
} else {
  ARGS <- list(max.edge.interior = 5, nres = 4, k = 2250, n.neighbours = 15, knots_squared = 20^2)
}

## Load the model fitting and prediction functions
invisible(
  mapply(source, paste0("./src/MODIS_modelling_fns/", PACKAGES, ".R"))
)

## Load user defined functions 
## These include functions to compute diagnostics, convert data frames from wide
## to long format in a convenient way for this study, and some other 
## miscellaneous functions to make the comparison study script easier to read.
source("./src/MODIS_diagnostic_and_misc_fns.R")

# ---- Load MODIS data  ----

data("MODIS_cloud_df") # MODIS dataframe stored in FRK (FRKTMB branch)
df <- MODIS_cloud_df; rm(MODIS_cloud_df)


# ---- Run comparison study and save results ----

number_of_runs_MAR <- 1
number_of_runs_block <- 1

## Missing at random
MAR_df_train_list <- MAR_df_test_list <- list()
MAR_times <- data.frame(matrix(ncol = length(PACKAGES) + 2, nrow = 0))
colnames(MAR_times) <- c(PACKAGES, "Sampling_scheme", "Run")
missing_fwk <- "MAR"
for (i in 1:number_of_runs_MAR) {
  seed <- i # Controls the training set
  source("./src/MODIS.R")
  
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
  seed <- i + 1 # NB: i + 1 gives an interesting block position when number_of_runs_block = 1
  source("./src/MODIS.R")
  
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
times <- times %>% gather(Method, time, PACKAGES)
times$Run <- as.integer(times$Run)
times$time <- as.numeric(times$time)
times$time <- times$time/60 # convert to minutes
write.csv(times,
          file = "./intermediates/times.csv",
          row.names = FALSE)



## ---- Read previously saved data ----

blocks <- read.csv("./intermediates/MODIS_blocks.csv")
all_df_train <- read.csv("./intermediates/MODIS_all_df_train.csv")
all_df_test <- read.csv(file = "./intermediates/MODIS_all_df_test.csv")
times <- read.csv(file = "./intermediates/times.csv")



## ---- Training and testing data visualisation ----

common_layers <- ggplot() + theme_bw() + coord_fixed() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + 
  labs(x = bquote(s[1]), y = bquote(s[2]))

## Plot the unthresholded version of the data
g_original_data <- common_layers +
  geom_raster(data = df, aes(x, y, fill = z_unthresholded)) +
  scale_fill_gradient(low = "black", high = "white") +
  labs(fill = "Radiance")

## Plot the thresholded data (which we use for our analyses)
discrete_cloud_scale <- scale_fill_gradient(low = "black", high = "white",
                                            breaks = c(0, 1),
                                            guide = "legend",
                                            labels = c("No Cloud", "Cloud"),
                                            name = "")

discrete_cloud_theme <- theme(legend.key=element_rect(fill = "gray", colour = "gray", size = 2))


training_data_background <- theme(
  panel.background = element_rect(fill = "#BFD5E3", colour = "#6D9EC1"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

## Full data set
g_thresholded_data <- common_layers +
  geom_raster(data = df, aes(x, y, fill = z)) +
  discrete_cloud_scale + discrete_cloud_theme


# ## Plot the training data.
# ## Just plot a single training set at a time; even if we use many for the final
# ## comparison, we will not include a plot of all sets, as it would take up 
# ## too much space.
# plot_training_testing_data <- function(df, scheme) {
#   common_layers +
#     geom_raster(data = df %>% 
#                   filter(Sampling_scheme == scheme, Run == 1), 
#                 aes(x, y, fill = z)) +
#     discrete_cloud_scale + discrete_cloud_theme + training_data_background
# }
# training_data_MAR_plot <- plot_training_testing_data(all_df_train, "MAR")
# test_data_MAR_plot <- plot_training_testing_data(all_df_test, "MAR")
# training_data_block_plot <- plot_training_testing_data(all_df_train, "block")
# suppressMessages(test_data_block_plot <- plot_training_testing_data(all_df_test, "block") +
#   lims(x = range(all_df_train$x), y = range(all_df_train$y)))
# 
# ## Now make common edits to the plots
# plot_list_data <- list(training_data_MAR_plot, test_data_MAR_plot, 
#                        training_data_block_plot, test_data_block_plot)
# 
# plot_list_data <- lapply(plot_list_data, function (gg) gg + theme(axis.text = element_text(size = 16),
#                                                                   axis.title = element_text(size = 19), 
#                                                                   legend.text = element_text(size = 16), 
#                                                                   plot.title = element_text(hjust = 0.5, size = 19)))
# 
# ggsave(
#   ggarrange(
#     plotlist = plot_list_data,
#     nrow = 2, ncol = 2, common.legend = TRUE, legend = "right"
#   ),
#   filename = "MODIS_data.png", device = "png", width = 10, height = 5.5,
#   path = "./img/"
# )

## With facetting:
tmp1 <- all_df_train %>% mutate(set = "training")
tmp2 <- all_df_test %>% mutate(set = "test", time = NULL)
tmp <- rbind(tmp1, tmp2) %>%
  mutate(facet_var = paste(Sampling_scheme, set)) %>%
  mutate(facet_var = factor(facet_var, 
                            levels = c("MAR training", "MAR test", "block training", "block test"),
                            labels = c("MR: training set", "MR: test set", "MB: training set", "MB: test set"))) %>%
  filter(Run == 1)

ggsave(
  common_layers +
    geom_raster(data = tmp, aes(x = x, y = y, fill = z))  +
    facet_wrap(~facet_var, nrow = 1) +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 13), 
          legend.text = element_text(size = 10), 
          strip.text = element_text(size = 12)) +
    discrete_cloud_scale + discrete_cloud_theme + training_data_background,
  filename = "MODIS_data.png", device = "png", width = 10, height = 2.2,
  path = "./img/"
)

# ## 2x2 facet grid:
# tmp1 <- all_df_train %>% mutate(set = "training")
# tmp2 <- all_df_test %>% mutate(set = "test", time = NULL)
# tmp <- rbind(tmp1, tmp2) 
# ggsave(
#   common_layers +
#     geom_raster(data = tmp, aes(x = x, y = y, fill = z))  +
#     facet_grid(Sampling_scheme ~ set) +
#     discrete_cloud_scale + discrete_cloud_theme + training_data_background,
#   filename = "MODIS_data.png", device = "png", width = 13, height = 2.2,
#   path = "./img/"
# )


# ---- Compute diagnostics (Brier score and AUC) ----

## Convert all_df_test to long format, and add run times
suppressMessages(
  all_df_test <- all_df_test %>% 
    gather(Method, pred, paste0("pred_", PACKAGES)) %>%
    mutate(Method = str_remove(Method, "pred_")) %>% 
    left_join(times)
)

## Compute diagnostics, split by run, method, and sampling scheme
suppressMessages(
  diagnostics <-  all_df_test %>%
    group_by(Sampling_scheme, Run, Method) %>%
    summarise(Brier = BrierScore(z, pred),
              AUC = AUC(z, pred), 
              time = time[1])
)

## Compute average diagnostics, split by method and sampling scheme, but averaged over all runs
suppressMessages(
  diagnostics <- diagnostics %>%  
    group_by(Sampling_scheme, Method) %>%
    summarise(
      Brier_mean = mean(Brier),
      # Brier_sd = sd(Brier),
      AUC_mean = mean(AUC), 
      # AUC_sd = sd(AUC), 
      time_mean = mean(time),
      # time_sd = sd(time)
    ) %>% 
    as.data.frame()
)


for (i in 1:ncol(diagnostics)) {
  if(is.numeric(diagnostics[, i])) 
    diagnostics[, i] <- round(diagnostics[, i], 3)
}


write.csv(diagnostics, row.names = FALSE, file = "./results/MODIS_diagnostics.csv")

stargazer::stargazer(
  diagnostics, summary = FALSE, rownames = FALSE,  
  out = "./results/MODIS_diagnostics.tex", 
  title = "Diagnostic results for the MODIS comparison study. Best performers for a given diagnostics are boldfaced", 
  label = "tab:MODIS_diagnostics")


# ---- Prediction maps ----

## Convert all_df_train to long format
all_df_train <- all_df_train %>% 
  gather(Method, pred, paste0("pred_", PACKAGES)) %>%
  mutate(Method = str_remove(Method, "pred_")) 

## Construct plotting data frame. Note that the location order does not
## matter for plotting, so we can just rbind all_df_train and all_df_test.
df_plot <- all_df_test %>%
  mutate(time = NULL) %>% 
  rbind(all_df_train)

## Add "Data" to Method: This represent the true, withheld data. 
N_test_locations <- nrow(all_df_test) / length(PACKAGES)
tmp <- all_df_test[1:N_test_locations, ]
tmp$Method <- "Data"
tmp$pred <- tmp$z
all_df_test <- rbind(all_df_test, tmp)

## Re-order the levels of Method to change the order of the panels
# df_plot$Method <- factor(df_plot$Method, levels = PACKAGES)

plot_predictions <- function(df_plot, scheme) {
  df_plot %>% 
    filter(Sampling_scheme == scheme, Run == 1) %>%
    ggplot() +
    geom_raster(aes(x, y, fill = pred)) +
    # facet_grid( ~ Method) +
    facet_wrap( ~ Method, nrow = 2) +
    scale_fill_gradient2(midpoint = 0.5, low = "black", mid = "gray", high = "white") +
    labs(fill = "") +
    theme_bw() + coord_fixed() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 13), 
          legend.text = element_text(size = 10), 
          strip.text = element_text(size = 12))  
}



ggsave(
  plot_predictions(all_df_test, "MAR") + training_data_background + 
    labs(x = bquote(s[1]), y = bquote(s[2])),
  filename = "MODIS_MAR_predictions.png", device = "png", 
  # width = 12, height = 2.35,
  width = 10, height = 5.1,
  path = "./img/"
)



## Focus in on the block region:
blocks <- read.csv("./intermediates/MODIS_blocks.csv")
ggsave(
  df_plot %>% 
    filter(Sampling_scheme == "block", Run == 1, Method == "FRK") %>%
    mutate(Method = "Data", pred = z) %>%
    rbind(df_plot) %>% 
    filter(Sampling_scheme == "block", Run == 1) %>%
    ggplot() +
    geom_raster(aes(x, y, fill = pred)) +
    facet_grid( ~ Method) +
    scale_fill_gradient2(midpoint = 0.5, low = "black", mid = "gray", high = "white") +
    labs(fill = "") +
    scale_x_continuous(limits = c(blocks$xmin[1], blocks$xmax[1]), expand = c(0, 0)) +
    scale_y_continuous(limits = c(blocks$ymin[1], blocks$ymax[1]), expand = c(0, 0)) +
    theme_bw() + coord_fixed() + 
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 13), 
          legend.text = element_text(size = 10), 
          strip.text = element_text(size = 12)) + 
    labs(x = bquote(s[1]), y = bquote(s[2])),
  filename = "MODIS_block_predictions.png", device = "png", width = 12, height = 2.5,
  path = "./img/"
)   

# ---- ROC curves ----

compute_ROC_objects <- function(df, scheme, run = 1) {
  ROC_list <- list()
  for (method in unique(all_df_test$Method)) {
    ROC_list[[method]] <- all_df_test %>%
      filter(Method == method) %>% 
      filter(Sampling_scheme == scheme, Run == run) %>% 
      pROC::roc(z, pred)
  }
  return(ROC_list)
}

suppressMessages({
  block_ROC_list <- compute_ROC_objects(all_df_test, scheme = "block")
  MAR_ROC_list <- compute_ROC_objects(all_df_test, scheme = "MAR")
})

FRK_idx <- which(PACKAGES == "FRK")

plot_block_ROC <- pROC::ggroc(
  block_ROC_list[c(PACKAGES[-FRK_idx], "FRK")], 
  legacy.axes = TRUE, size = 0.3) +
  scale_colour_manual(values = 1:5) +
  xlab("False positive rate") + ylab("True positive rate") +
  labs(colour = "Method") +
  geom_abline(slope = 1, intercept = 0, color = "darkgrey", linetype = "dashed") +
  theme_bw() + coord_fixed()

plot_MAR_ROC <- pROC::ggroc(
  MAR_ROC_list[c(PACKAGES[-FRK_idx], "FRK")],
  legacy.axes = TRUE, size = 0.3) +
  scale_colour_manual(values = 1:5) +
  xlab("False positive rate") + ylab("True positive rate") +
  labs(colour = "Method") +
  geom_abline(slope = 1, intercept = 0, color = "darkgrey", linetype = "dashed") +
  theme_bw() + coord_fixed()


ggsave(
  ggpubr::ggarrange(plot_MAR_ROC, plot_block_ROC, align = "h",
                    common.legend = TRUE, legend = "right"),
  filename = "MODIS_ROC.png", device = "png", width = 6, height = 2.6, # plot together
  path = "./img/"
)

