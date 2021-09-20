library("FRK")
library("plyr")     # round_any()
library("dplyr")    
library("ggplot2")
library("sp")
library("maptools") # readShapePoly()
library("raster")   # bind()
library("ggpubr")
library("ggmap")    # Stamenmap

source("./scripts/Utility_fns.R")

## Use low-rank versions of the models to establish that the code works? 
quick <- check_quick()
nres <- if (quick) 2 else 3


## Some helper functions for plotting
source("./scripts/Plotting_helpers/Sydney_plotting_fns.R")
source("./scripts/Plotting_helpers/Plotting_helpers.R")
source("./scripts/Utility_fns.R")

# ---- Data preprocessing ----

source("./scripts/Sydney_prep.R")

# ---- Sydney Stamen map ----

## map background to show Sydney
## NB: get_stamenmap() does NOT require a google API
Sydney_bbox <- c(left = 150.72, bottom = -34.2, right = 151.32, top = -33.65)
Sydney_map  <- get_stamenmap(bbox = Sydney_bbox, 
                             maptype = "toner-background", 
                             color = "bw")
save(Sydney_map, file = "./data/Sydney_map.RData") 

## Create map layer to place under all plots
Sydney_map <- ggmap(Sydney_map)


# ---- Conduct the analysis ---- 

## The argument fitting controls whether the training data comprises SA2 regions 
## only (fitting = "SA2s"), SA1 regions only (fitting = "SA1s"; this is useful 
## to obtain a good estimate of the fine-scale variance parameter), or a mixture 
## of SA2s and SA1s (fitting = "mixed"; this is the approach used in the paper).  
Sydney_analysis <- function(fitting = "mixed") {
  
  # ---- Construct the training data ----
  
  ## Data for model fitting
  poly_fit <- construct_training_data(fitting = fitting, SA1s = SA1s, SA2s = SA2s)
  
  ## Remove SA regions that have no families of interest
  poly_fit <- subset(poly_fit, number_of_families > 0)
  
  # ---- Plotting training data ----
  
  ## gray background to show SA2s with no families of interest
  SA2_bg <-  geom_polygon(
    data = SpatialPolygonsDataFrame_to_df(SA2_NSW_sub_withk0), 
    aes(lon, lat, group = SA2_MAIN11), fill = "light gray"
  )
  
  training_data_plots <- plot_spatial_or_ST(
    poly_fit, 
    c("number_of_families", "Proportion_poverty"),
    map_layer = Sydney_map + SA2_bg, colour = "black", size = 0.1
  )
  
  ## Adjust the legend breaks and legend labels:
  breaks <- list(
    number_of_families = c(2000, 6000), 
    Proportion_poverty = if (fitting == "SA2s") c(0.1, 0.2, 0.3) else c(0.25, 0.75)
  )
  
  fill_label = list(
    number_of_families = expression(paste("Number of families, ", bold(k)[Z])),
    Proportion_poverty  = "Observed proportion\nof families in poverty"
  )
  
  ## axis labels (used for all plots)
  lab1 <- xlab("lon (deg)")
  lab2 <- ylab("lat (deg)")
  
  training_data_plots <- lapply(
    names(training_data_plots), function(i) {
      gg <- training_data_plots[[i]] 
      gg <- gg %>% change_legend_breaks(breaks[[i]]) 
      gg <- gg + labs(fill = fill_label[[i]]) + lab1 + lab2
      gg <- gg %>% change_font_size_and_axis
      return(gg)
    }
  )
  
  ggsave( 
    ggarrange(plotlist = training_data_plots, align = "hv", nrow = 1, legend = "top"),
    filename = paste0("4_3_Sydney_training_data_", fitting, ".png"), 
    device = "png", path = "./results/", width = 9, height = 4.1
  )
  
  # ---- FRK ----
  
  BAUs <- SA1s
  
  ## Assign the size parameter of each BAU
  BAUs$k_BAU <- BAUs$number_of_families
  
  ## Remove overlapping columns from the BAUs and the data:
  ## (i.e., number_of_families_in_poverty, number_of_families, etc.).
  BAUs@data[, which(names(BAUs) %in% names(poly_fit))] <- NULL
  BAUs$fs <- 1 # homoscedastic fine-scale variation
  
  ## For binomial or negative-binomial data, the known constant parameter k must
  ## be provided for each observation. In this case, it is simply the total number
  ## of families of interest
  poly_fit$k_Z <- poly_fit$number_of_families
  
  if (fitting == "SA2s") {
    ## A reliable estimate of the fine-scale variance can be obtained by first 
    ## fitting the model using SA1 data (in practice, one may use past census data).
    known_sigma2fs <- 0.2127 
  } else {
    ## Estimate fine-scale variance using TMB (known_sigma2fs = NULL)
    known_sigma2fs <- NULL
  }
  
  ## Construct and fit the SRE object
  S <- FRK(f = number_of_families_in_poverty ~ 1,
           nres = nres,
           data = list(poly_fit), BAUs = BAUs, 
           response = "binomial", link = "logit", 
           known_sigma2fs = known_sigma2fs) 
  
  # ---- Prediction over the SA1s ----
  
  RNGversion("3.6.0"); set.seed(1)
  pred <- predict(S, type = c("mean", "response"))
  
  # ---- Model validation using SA1 level data ----
  
  ## Coverage
  val_id     <- which(SA1_NSW_sub$number_of_families > 10)
  true_value <- SA1_NSW_sub[val_id, ]$number_of_families_in_poverty
  lower      <- pred$newdata@data[val_id, "Z_percentile_5"]
  upper      <- pred$newdata@data[val_id, "Z_percentile_95"]
  coverage   <- mean((lower <= true_value) & (true_value <= upper))
  
  diagnostics <- data.frame(coverage = coverage)
  save_path   <- paste0("./results/4_3_Sydney_SA1_coverage_", fitting)
  
  write.csv(diagnostics, 
            paste0(save_path, ".csv"), 
            row.names = FALSE)
  
  save_html_table(
    diagnostics,
    file = paste0(save_path, ".html"), 
    caption = "Sydney spatial change-of-support"
  )
  
  
  
  # ---- Plotting SA1 predictions ----
  
  plots <- plot(
    S, pred$newdata, 
    map_layer = Sydney_map + SA2_bg,             # optional layer to put below the plotting geom
    colour = "black", size = 0.025, alpha = 0.9  # optional arguments to plotting geom via ...
  ) 
  
  ## Change axis labels and font size
  plots <- lapply(plots, function(gg) gg + lab1 + lab2) 
  plots <- lapply(plots, change_font_size_and_axis)
  
  ## Simple legends (useful for presentations)
  if (fitting == "SA2s") plots <- simplify_legend_label(plots)
  
  ## Shift the legend title so it doesn't run into the legend box
  plots$interval90_prob <- plots$interval90_prob + theme(legend.title.align = 1.5)
  
  ggsave( 
    ggarrange(plots$p_prob, plots$interval90_prob, nrow = 1, legend = "top"),
    filename = paste0("4_3_Sydney_SA1_predictions_", fitting, ".png"), 
    device = "png", width = 9.5, height = 4.4, path = "./results/"
  )
  
  # ---- SA3 predictions ----
  
  RNGversion("3.6.0"); set.seed(1)
  pred <- predict(S, newdata = SA3_NSW_sub)
  plots <- plot(S, pred$newdata, colour = "black", map_layer = Sydney_map, alpha = 0.9)
  
  ## Change the breaks, labels, and font size
  breaks <- list(
    p_prob          = c(0.1, 0.16, 0.22), 
    interval90_prob = c(0.006, 0.009, 0.012), 
    p_mu            = c(2000, 6000, 10000), 
    interval90_mu   = c(150, 225, 300)
  )
  
  plotnames <- names(plots)
  plots <- lapply(
    names(plots), function(i) {
      gg <- plots[[i]] 
      gg <- gg %>% change_legend_breaks(breaks[[i]]) 
      gg <- gg + lab1 + lab2
      gg <- gg %>% change_font_size_and_axis
      return(gg)
    }
  )
  names(plots) <- plotnames
  
  ## Simple legends (useful for presentations)
  if (fitting == "SA2s") plots <- simplify_legend_label(plots)
  
  ## Shift the legend title so it doesn't run into the legend box
  plots$interval90_prob <- plots$interval90_prob + theme(legend.title.align = 1.8)
  plots$interval90_mu <- plots$interval90_mu + theme(legend.title.align = 1.5)
  
  ggsave( 
    ggarrange(plots$p_prob, plots$interval90_prob, nrow = 1, legend = "top"),
    filename = paste0("4_3_Sydney_SA3_predictions_probability_", fitting, ".png"), 
    device = "png", width = 9.5, height = 4.4, path = "./results/"
  )
}

cat("Conducting Sydney analysis: Using a mixture of SA1 and SA2 regions as training data")
Sydney_analysis(fitting = "mixed")

## This is only used for presentations, and not at all in the paper:
# cat("Conducting Sydney analysis: Using only SA2 regions as training data")
# Sydney_analysis(fitting = "SA2s")