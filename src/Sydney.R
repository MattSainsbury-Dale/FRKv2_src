library("FRK")
library("plyr")     # round_any()
library("dplyr")    
library("ggplot2")
library("sp")
library("maptools") # readShapePoly()
library("raster") # bind()
library("ggpubr")
library("ggmap") # Stamenmap

defaultW <- getOption("warn") 
options(warn = -1)

# ---- Data preprocessing ----

source("./src/Sydney_prep.R")

# ---- Sydney Stamen map ----

## map background to show Sydney
## NB: get_stamenmap() does NOT require a google API.
## NB: To revert to original format, just remove Sydney from map_layer (in all plots), 
## and remove expand = c(0, 0) in change_font_size_and_axis(). 
Sydney_bbox = c(left = 150.72, bottom = -34.2, right = 151.32, top = -33.65)
Sydney_map <- get_stamenmap(bbox = Sydney_bbox, 
                            maptype = "toner-background", color = "bw")
save(Sydney_map, file = "./data/Sydney_map.RData") 

## Create map layer to place under all plots
Sydney_map <- ggmap(Sydney_map)


## Conducts the analysis. The argument fitting controls whether the training data
## consists of SA2 regions only (fitting = "SA2s"), SA1 regions only (fitting = 
## "SA1s"; useful to get a good estimate of the fine-scale variance parameter), 
## or a mixture of SA2s and SA1s (fitting = "mixed"; this is the approach used in 
## the paper).  
Sydney_analysis <- function(fitting = "mixed") {
  
  # ---- Construct the training data ----
  
  ## Data for model fitting
  poly_fit <- construct_training_data(fitting = fitting, SA1s = SA1s, SA2s = SA2s)
  
  ## Remove SA regions that have no families of interest
  poly_fit <- subset(poly_fit, Total_families_of_interest > 0)

  # ---- Plotting training ----
  
  ## ggplot annotations
  lab1 <- xlab("lon (deg)")
  lab2 <- ylab("lat (deg)")
  
  ## gray background to show SA2s with no families of interest
  SA2_bg <-  geom_polygon(data = SpatialPolygonsDataFrame_to_df(SA2_NSW_sub_withk0),
                          aes(lon, lat, group = SA2_MAIN11),
                          fill = "light gray")
  
  column_names <- c("Total_families_of_interest", "Proportion_poverty")
  
  training_data_plots <- plot_spatial_or_ST(
    poly_fit, column_names, map_layer = Sydney_map + SA2_bg, 
    colour = "black", size = 0.1
  )
  
  fill_label = c(
    Total_families_of_interest = expression(paste("Number of families, ", bold(k)[Z])),
    Proportion_poverty  = "Observed proportion\nof families in poverty"
  )
  
  
  if (fitting == "SA2s") {
    breaks <- list(c(2000, 6000), c(0.1, 0.2, 0.3)) 
  } else {
    breaks <- list(c(2000, 6000), c(0.25, 0.75))
  }
  
  training_data_plots <- lapply(1:length(training_data_plots), 
                                function(i) training_data_plots[[i]] + 
                                  labs(fill = fill_label[i]) + lab1 + lab2 + 
                                  scale_fill_distiller(palette = "Spectral", 
                                                       breaks = breaks[[i]]))
  
  names(training_data_plots) <- column_names
  
  change_font_size_and_axis <- function(gg) {
    gg + theme(axis.text = element_text(size = 11),
               axis.title = element_text(size = 14), 
               legend.text = element_text(size = 11),
               legend.title = element_text(size = 14)) + 
      scale_x_continuous(breaks = c(150.8, 151.0, 151.2), expand = c(0, 0)) + 
      scale_y_continuous(breaks = c(-33.7, -33.9, -34.1), expand = c(0, 0))
  }
  
  training_data_plots <- lapply(training_data_plots, change_font_size_and_axis)
  
  ggsave( 
    ggarrange(plotlist = training_data_plots, align = "hv", nrow = 1, legend = "top"),
    filename = paste0("Sydney_training_data_", fitting, ".png"), 
    device = "png", width = 9, height = 4.1,
    path = "./img/"
  )
  
  
  # ---- FRK ----
  
  BAUs <- SA1s
  
  ## Assign the size parameter of each BAU
  BAUs$k_BAU <- BAUs$Total_families_of_interest
  
  ## Remove overlapping columns from the BAUs and the data:
  ## (i.e., total_poverty_count, Total_families_of_interest, etc.).
  BAUs@data[, which(names(BAUs) %in% names(poly_fit))] <- NULL
  BAUs$fs <- 1 # homoscedastic fine-scale variation
  
  ## For binomial or negative-binomial data, the known constant parameter k must
  ## be provided for each observation. In this case, it is simply the total number
  ## of families of interest
  poly_fit$k_Z <- poly_fit$Total_families_of_interest
  
  
  if (fitting == "SA2s") {
    ## A reliable estimate of the fine-scale variance can be obtained by first 
    ## fitting the model using SA1 data (in practice, one may use past census data).
    known_sigma2fs <- 0.2127 
  } else {
    known_sigma2fs <- NULL
  }
  
  ## Construct and fit the SRE object
  S <- FRK(
    f = total_poverty_count ~ 1, data = list(poly_fit), BAUs = BAUs, 
    response = "binomial", link = "logit", 
    normalise_wts = FALSE,     # sum (not average) the mean-process over the SA1s
    known_sigma2fs = known_sigma2fs
  ) 
  
  # ---- Prediction over the SA1s ----
  
  RNGversion("3.6.0")
  set.seed(1)
  pred <- predict(S, type = c("mean", "response"))
  
  # ---- Model validation using SA1 level data ----
  
  ## Coverage
  val_id     <- which(SA1_NSW_sub$Total_families_of_interest > 10)
  true_value <- SA1_NSW_sub[val_id, ]$total_poverty_count
  lower      <- pred$newdata@data[val_id, "Z_percentile_5"]
  upper      <- pred$newdata@data[val_id, "Z_percentile_95"]
  coverage   <- mean((lower <= true_value) & (true_value <= upper))
  
  write.csv(data.frame(coverage = coverage), 
            paste0("./results/Sydney_SA1_coverage_", fitting, ".csv"), 
            row.names = FALSE)
  
  
  # ---- Plotting SA1 predictions ----
  
  ## We can pass in alpha via ..., however there is no way to set alpha for colour
  ## and fill to be different. One workaround is to provide an rgb specification 
  ## for colour, with a string of the form "#RRGGBB" where each of the pairs RR, 
  ## GG, BB consists of two hexadecimal digits giving a value in the range 00 to FF. 
  ## You can optionally make the colour transparent by using the form "#RRGGBBAA", 
  ## where AA denotes the alpha (opacity) level.
  ## See https://hashnode.com/post/understanding-rrggbbaa-color-notation-cisvdr52x088fwt53h1drf6m2
  ## for a useful blog post on the rgba notation.
  
  plots <- plot(
    S, pred$newdata, 
    map_layer = Sydney_map + SA2_bg,  # optional layer to put below the plotting geom
    colour = "black", size = 0.025    # optional arguments to plotting geom via ...
  ) 
  
  plots <- lapply(plots, function(gg) gg + lab1 + lab2)
  plots <- lapply(plots, change_font_size_and_axis)
  
  ## We use SA2s for presentations: Edit the legend labels so that they are as 
  ## simple as possible.
  if (fitting == "SA2s") {
    plots$p_prob <- plots$p_prob + 
      labs(fill = "Predicted proportion\nof families in poverty")
    plots$interval90_prob <- plots$interval90_prob + 
      # labs(fill = "90% prediction-interval width\nfor proportion of families in poverty")
      labs(fill = "90% prediction-interval \nwidth for proportion of\nfamilies in poverty")
  }
  
  ggsave( 
    ggarrange(plots$p_prob, plots$interval90_prob, nrow = 1, legend = "top"),
    filename = paste0("Sydney_SA1_predictions_", fitting, ".png"), 
    device = "png", width = 9.5, height = 4.4, path = "./img/"
  )
  
  
  # ---- SA3 predictions ----
  
  RNGversion("3.6.0"); set.seed(1)
  pred <- predict(S, newdata = SA3_NSW_sub)
  plots <- plot(S, pred$newdata, colour = "black", map_layer = Sydney_map)
  plots <- lapply(plots, function(gg) gg + lab1 + lab2)
  plots <- lapply(plots, change_font_size_and_axis)
  plots$p_mu <-  plots$p_mu + scale_fill_distiller(palette="Spectral", breaks = c(2000, 6000, 10000))
  plots$interval90_mu <-  plots$interval90_mu + scale_fill_distiller(palette="BrBG", breaks = c(150, 225, 300)) 
  
  ## We use SA2s for presentations: Edit the legend labels so that they are as 
  ## simple as possible.
  if (fitting == "SA2s") {
    plots$p_mu <- plots$p_mu + 
      labs(fill = "Predicted number\nof families in poverty")
    plots$interval90_mu <- plots$interval90_mu + 
      labs(fill = "90% prediction-interval width\nfor number of families in poverty")
  }
  
  
  ggsave( 
    ggarrange(plots$p_mu, plots$interval90_mu, nrow = 1, legend = "top"),
    filename = paste0("Sydney_SA3_predictions_", fitting, ".png"), 
    device = "png", width = 9.5, height = 4.4, path = "./img/"
  )
}


cat("Conducting Sydney analysis: Using a mixture of SA1 and SA2 regions as training data")
Sydney_analysis(fitting = "mixed")

cat("Conducting Sydney analysis: Using only SA2 regions as training data")
Sydney_analysis(fitting = "SA2s")

options(warn = defaultW)



