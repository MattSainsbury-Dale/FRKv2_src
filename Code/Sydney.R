suppressMessages({
library("FRK")
library("plyr")     # round_any()
library("dplyr")    
library("ggplot2")
library("sp")
library("sf")       # st_read() and st_is_empty()
library("raster")   # bind()
library("ggpubr")
library("ggmap")    # Stamenmap

source("Code/Utility_fns.R")
})

## Use very-low-dimensional representations of the models to establish that the code works? 
quick <- check_quick()
nres <- if (quick) 2 else 3

# ---- Data preprocessing ----

## Load the census data
census_SA1_df <- read.csv("data/Sydney_SA1_data.csv")
census_SA2_df <- read.csv("data/Sydney_SA2_data.csv")

## Poverty lines as defined in Appendix E: "Sydney poverty lines"
poverty_lines <- c("Couple_family_with_no_children" = round_any(594.6, 200), 
                   "Couple_family_with_children" = round_any(835.3, 200), 
                   "One_parent_family" = round_any(691.04, 200))

## Select the income brackets common to all family groups as defined above
low_income <- c("Negative_Nil_income_", "X1_199_", "X200_299_", "X300_399_", "X400_599_")

## Select the family names we are interested in.
family_names <- names(poverty_lines)

## The following function performs data pre-processing on the SA1 or SA2 data frames.
## At each SA, it first determines the number of families below the poverty line 
## within each family type, then computes the total number of families in poverty
## (i.e., the number of "successes", Z) and the number of families in total 
## (i.e., the size parameter, k). For the sake of data visualisation, we also 
## compute the empirical proportion of families in poverty.
censusDataPreprocess <- function(df, SA) {
  
  df <- df %>% 
    mutate(
      poverty_count_Couple_family_with_no_children = rowSums(.[, paste0(low_income, family_names[1])]),
      poverty_count_Couple_family_with_children = rowSums(.[, paste0(low_income, family_names[2])]), 
      poverty_count_One_parent_family = rowSums(.[, paste0(c(low_income, "X600_799_"), family_names[3])])
      
    ) %>% 
    mutate(
      number_of_families_in_poverty = poverty_count_Couple_family_with_no_children + 
        poverty_count_Couple_family_with_children + 
        poverty_count_One_parent_family
    ) %>% 
    mutate(
      number_of_families = 
        Total_Couple_family_with_no_children + 
        Total_Couple_family_with_children + 
        Total_One_parent_family
    )  %>% 
    mutate(
      Proportion_poverty = number_of_families_in_poverty / number_of_families
    ) %>% 
    dplyr::select(region_id, number_of_families_in_poverty, number_of_families, Proportion_poverty)
  
  return(df)
}

## Compute the total number of families, number of families in poverty, 
## and the proportion of families in poverty at SA1 and SA2 level. 
census_SA1_df <- censusDataPreprocess(census_SA1_df, "SA1")
census_SA2_df <- censusDataPreprocess(census_SA2_df, "SA2")


## SA Shapefiles and pre-processing

## Since SA3s are aggregations of SA2s, which in turn are
## aggregations of SA1s, we define the domain of interest first in terms of 
## SA3s, and then subset the SA2s and SA1s accordingly. This is to ensure that 
## the prediction maps look consistent, in particular, the SA3 predictions don't 
## have a big chunk of SA1s missing.  This will ensure that the SA1s and SA2s 
## used for model fitting aggregate to whole SA3s. 

## Helper function for reading spatial polygons saved as a .shp file
read_shape_poly <- function(path) {
  
  # polys <- maptools::readShapePoly(path, delete_null_obj = TRUE)
  
  polys <- sf::st_read(path, quiet = TRUE) %>%
    filter(!sf::st_is_empty(.)) %>%
    as("Spatial")
  
  crs(polys) <- NA
  
  return(polys)
}

## Deal with the SA3s first so that we can subset the SA1s based on the SA3s. 
SA3 <- read_shape_poly("data/Sydney_shapefiles/SA3/SA3_2011_AUST.shp")
coordnames(SA3) <- c("lon", "lat")

## Now focus on a subset of SA3s in NSW (around Sydney)
cc          <- coordinates(SA3)
idx         <- cc[, "lon"] > 150.77 & cc[, "lon"] < 151.5 & cc[, "lat"] > -34.1 & cc[, "lat"] < -33.73
SA3_NSW_sub <- SA3[idx,]
# number of SA3s used for fitting
cat("Number of SA3s:", length(SA3_NSW_sub), "\n")

## Deal with the SA2s first so that we can subset the SA1s based on the SA2s 
SA2 <- read_shape_poly("data/Sydney_shapefiles/SA2/SA2_2011_AUST.shp")
coordnames(SA2) <- c("lon", "lat")

## Add the SA2 census data 
census_SA2_df$SA2_MAIN11 <- census_SA2_df$region_id
SA2@data <- merge(SA2@data, census_SA2_df, by = "SA2_MAIN11")

## Retain only those SA2s which are in the remaining SA3s
SA2_NSW_sub <- subset(SA2, SA3_CODE11 %in% SA3_NSW_sub$SA3_CODE11) 

## Remove any SA2s which had 0 total families of interest
SA2_NSW_sub_withk0 <- SA2_NSW_sub  # keep a backup with all of the SA2s
idx <- which(SA2_NSW_sub@data$number_of_families == 0)
cat("SA2s with no families of interest: ", 
    paste(SA2_NSW_sub@data$SA2_NAME11[idx], collapse=", "), "\n")
SA2_NSW_sub <- subset(SA2_NSW_sub, number_of_families > 0)
cat("Removed SA2s with no families of interest.\n")
cat("Number of SA2s:", length(SA2_NSW_sub), "\n")

## Now deal with the SA1s
SA1 <- read_shape_poly("data/Sydney_shapefiles/SA1/SA1_2011_AUST.shp")
coordnames(SA1) <- c("lon", "lat")

## Retain only those SA1s which are in the remaining SA2s.
SA1_NSW_sub <- subset(SA1, SA2_MAIN11 %in% SA2_NSW_sub$SA2_MAIN11) 
cat("Number of SA1s:", length(SA1_NSW_sub), "\n")

## Add the SA1 census data
census_SA1_df$SA1_7DIG11 <- census_SA1_df$region_id
SA1_NSW_sub@data <- merge(SA1_NSW_sub@data, census_SA1_df,
                          by = "SA1_7DIG11")

## For code simplicty, remove _NSW_sub henceforth
SA1s <- SA1_NSW_sub
SA2s <- SA2_NSW_sub

construct_training_data <- function(fitting = "mixed", SA1s, SA2s) {
  if (fitting == "SA2s") {
    return(SA2s)
  } else if (fitting == "SA1s") {
    return(SA1s)
  } else if (fitting == "mixed") {
    
    ## We include some SA1s in the training set 
    zdf  <- SA2s 
    
    ## Select some data supports to be replaced
    RNGversion("3.6.0")
    set.seed(123)
    rmidx <- sample(1:length(zdf), ceiling(length(zdf)/10), replace = FALSE)
    
    ## Convert SA1s to points, otherwise SA1s in neighbouring data supports are 
    ## also selected because they share lines, which is classified as overlapping 
    ## by over()
    SA1s_as_points <- SpatialPoints(as(SA1s, "SpatialPolygons"))
    crs(SA1s_as_points) <- crs(zdf)
    ## Select SA1s corresponding to the selected data supports
    observed_BAU_idx <- which(!is.na(over(SA1s_as_points, zdf[rmidx, 1])))
    
    ## Do not add SA1 regions that have less than 0 families of interest (this would
    ## cause complications in model fitting)
    observed_BAU_idx <- intersect(observed_BAU_idx, 
                                  which(SA1s$number_of_families > 0))
    
    ## Combine the observed SA1s with the SA2s
    zdf <- bind(SA1s[observed_BAU_idx, ], zdf[-rmidx, ])
    
    return(zdf)
  }
}


# ---- Sydney Stamen map ----

## Note that get_stamenmap() does not currently require a google API, but I 
## chose to save the map object in case this changes in the future
# Sydney_bbox <- c(left = 150.72, bottom = -34.2, right = 151.32, top = -33.65)
# Sydney_map  <- get_stamenmap(bbox = Sydney_bbox, 
#                              maptype = "toner-background", 
#                              color = "bw")
# save(Sydney_map, file = "data/Sydney_map.RData") 

## Create map layer to place under all plots
load("data/Sydney_map.RData") 
Sydney_map <- ggmap(Sydney_map)


# ---- Conduct the analysis ---- 

## The argument `fitting` controls whether the training data comprise SA2 regions 
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
  
  suppressMessages({
  training_data_plots <- plot_spatial_or_ST(
    poly_fit, 
    c("number_of_families", "Proportion_poverty"),
    map_layer = Sydney_map + SA2_bg, colour = "black", size = 0.1
  )
  })
  
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
  
  suppressMessages({
  training_data_plots <- lapply(
    names(training_data_plots), function(i) {
      gg <- training_data_plots[[i]] 
      gg <- gg %>% change_legend_breaks(breaks[[i]]) 
      gg <- gg + labs(fill = fill_label[[i]]) + lab1 + lab2
      return(gg)
    }
  )
  })
  
  figure <- ggarrange(plotlist = training_data_plots, align = "hv", nrow = 1, legend = "top")
  
  ggsave( 
    ggarrange(plotlist = training_data_plots, align = "hv", nrow = 1, legend = "top"),
    filename = paste0("4_3_Sydney_training_data_", fitting, ".png"), 
    device = "png", path = "Figures/", width = 9, height = 4.1
  )
  
  # ---- FRK ----
  
  BAUs <- SA1s
  
  ## For binomial or negative-binomial data, the known constant parameter k must
  ## be provided for each observation. In this case, it is simply the total number
  ## of families of interest. However, where we declare it depends on whether 
  ## the observation supports are all contained within BAUs, or if some observation
  ## supports comprise multiple BAUs.
  
  if (fitting == "SA1s") {
    ## Estimate fine-scale variance using TMB (known_sigma2fs = NULL)
    known_sigma2fs <- NULL
    poly_fit$k_Z <- poly_fit$number_of_families
    
    ## Could also provide the size parameter like this: 
    # BAUs$k_BAU <- BAUs$number_of_families

  } else if (fitting == "SA2s") {
    ## A reliable estimate of the fine-scale variance can be obtained by first 
    ## fitting the model using SA1 data (in practice, one may use past census data).
    known_sigma2fs <- 0.2127 
    BAUs$k_BAU <- BAUs$number_of_families
  }  else if (fitting == "mixed") {
    ## Estimate fine-scale variance using TMB (known_sigma2fs = NULL)
    known_sigma2fs <- NULL
    BAUs$k_BAU <- BAUs$number_of_families
  }
  
  ## Remove overlapping columns from the BAUs and the data:
  ## (i.e., number_of_families_in_poverty, number_of_families, etc.).
  BAUs@data[, which(names(BAUs) %in% names(poly_fit))] <- NULL
  BAUs$fs <- 1 # homoscedastic fine-scale variation

  ## Construct and fit the SRE object
  S <- FRK(f = number_of_families_in_poverty ~ 1,
           nres = nres,
           data = list(poly_fit), BAUs = BAUs, 
           response = "binomial", link = "logit", 
           known_sigma2fs = known_sigma2fs, 
           # manually set these arguments to reduce console output:
           K_type = "precision", method = "TMB", est_error = FALSE) 
  
  # ---- Prediction over the SA1s ----
  
  RNGversion("3.6.0"); set.seed(1)
  pred <- predict(S, type = c("mean", "response"), k = BAUs$number_of_families)
  
  # ---- Model validation using SA1 level data ----
  
  ## Coverage
  val_id     <- which(SA1_NSW_sub$number_of_families > 10)
  true_value <- SA1_NSW_sub[val_id, ]$number_of_families_in_poverty
  lower      <- pred$newdata@data[val_id, "Z_percentile_5"]
  upper      <- pred$newdata@data[val_id, "Z_percentile_95"]
  coverage   <- mean((lower <= true_value) & (true_value <= upper))
  
  diagnostics <- data.frame(coverage = coverage)
  save_path   <- paste0("Figures/4_3_Sydney_SA1_coverage_", fitting)
  
  write.csv(diagnostics, 
            paste0(save_path, ".csv"), 
            row.names = FALSE)
  
  save_html_table(
    diagnostics,
    file = paste0(save_path, ".html"), 
    caption = "Sydney spatial change-of-support"
  )
  
  # ---- Plotting SA1 predictions ----
  
  suppressMessages({
  plots <- plot(
    S, pred$newdata, 
    map_layer = Sydney_map + SA2_bg,             # optional layer to put below the plotting geom
    colour = "black", size = 0.025, alpha = 0.9  # optional arguments to plotting geom via ...
  ) 
  
  ## Change axis labels and font size
  plots <- lapply(plots, function(gg) gg + lab1 + lab2) 
  
  ## Shift the legend title so it doesn't run into the legend box
  plots$interval90_prob <- plots$interval90_prob + theme(legend.title.align = 1.5)
  })
  
  suppressWarnings(
    figure <- ggarrange(plots$p_prob, plots$interval90_prob, nrow = 1, legend = "top")
    )
  
  ggsave( 
    figure,
    filename = paste0("4_3_Sydney_SA1_predictions_", fitting, ".png"), 
    device = "png", width = 9.5, height = 4.4, path = "Figures/"
  )
  
  # ---- SA3 predictions ----
  
  RNGversion("3.6.0"); set.seed(1)
  pred <- predict(S, newdata = SA3_NSW_sub)
  
  ## breaks, labels, and font size
  breaks <- list(
    p_prob          = c(0.1, 0.16, 0.22), 
    interval90_prob = c(0.006, 0.009, 0.012), 
    p_mu            = c(2000, 6000, 10000), 
    interval90_mu   = c(150, 225, 300)
  )
  
  suppressMessages({
  plots <- plot(S, pred$newdata, colour = "black", map_layer = Sydney_map, alpha = 0.9)
  plotnames <- names(plots)
  plots <- lapply(
    names(plots), function(i) {
      gg <- plots[[i]] 
      gg <- gg %>% change_legend_breaks(breaks[[i]]) 
      gg <- gg + lab1 + lab2
      return(gg)
    }
  )
  })
  names(plots) <- plotnames
  
  ## Shift the legend title so it doesn't run into the legend box
  plots$interval90_prob <- plots$interval90_prob + theme(legend.title.align = 1.8)
  plots$interval90_mu <- plots$interval90_mu + theme(legend.title.align = 1.5)

  suppressWarnings(
    figure <- ggarrange(plots$p_prob, plots$interval90_prob, nrow = 1, legend = "top")
  )
  
  ggsave( 
    figure,
    filename = paste0("4_3_Sydney_SA3_predictions_probability_", fitting, ".png"), 
    device = "png", width = 9.5, height = 4.4, path = "Figures/"
  )
}

cat("\nStarting FRK analysis: Using a mixture of SA1 and SA2 regions as training data.\n")
Sydney_analysis(fitting = "mixed")

## This is only used for presentations, and not at all in the paper:
# cat("\nStarting FRK analysis: Using only SA2 regions as training data.\n")
# Sydney_analysis(fitting = "SA2s")

## This serves a test for the case of binomial data with areal observation
## supports, where the observation supports and BAUs coincide:
# cat("\nStarting FRK analysis: Using only SA1 regions as training data.\n")
# Sydney_analysis(fitting = "SA1s")
