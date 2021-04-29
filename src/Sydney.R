library("FRK")
library("plyr")     # round_any()
library("dplyr")    
library("ggplot2")
library("sp")
library("maptools") # readShapePoly()
library("raster") # bind()
library("ggpubr")

defaultW <- getOption("warn") 
options(warn = -1)

# ---- Data preprocessing ----

source("./src/Sydney_prep.R")

# ---- Construct the training data ----

## Data for model fitting
poly_fit <- construct_training_data(fitting = "mixed", SA1s = SA1s, SA2s = SA2s)

## Remove SA regions that have no families of interest
poly_fit <- subset(poly_fit, Total_families_of_interest > 0)


# ---- Plotting training ----

## ggplot annotations
lab1 <- xlab("lon (deg)")
lab2 <- ylab("lat (deg)")

## gray background to show SA2s with no families of interest
SA2_bg <- geom_polygon(data = SpatialPolygonsDataFrame_to_df(SA2_NSW_sub_withk0),
                       aes(lon, lat, group = SA2_MAIN11),
                       fill = "light gray")

column_names <- c("Total_families_of_interest", "Proportion_poverty")
training_data_plots <- plot_spatial_or_ST(poly_fit, column_names, map_layer = SA2_bg, colour = "black", size = 0.1)

fill_label = c(
  Total_families_of_interest = "Number of families (k)",
  Proportion_poverty  = "Observed proportion\nof families in poverty"
)

breaks <- list(c(2000, 6000), c(0.25, 0.75))

suppressMessages(
  training_data_plots <- lapply(1:length(training_data_plots), 
                                function(i) training_data_plots[[i]] + 
                                  labs(fill = fill_label[i]) + lab1 + lab2 + 
                                  scale_fill_distiller(palette = "Spectral", breaks = breaks[[i]])
  )
)
names(training_data_plots) <- column_names

ggsave( 
  ggarrange(plotlist = training_data_plots, align = "hv", nrow = 1, legend = "top"),
  filename = "Sydney_training_data.png", device = "png", width = 9, height = 4.1,
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

## Construct and fit the SRE object
S <- FRK(
  f = total_poverty_count ~ 1, data = list(poly_fit), BAUs = BAUs, 
  response = "binomial", link = "logit", 
  normalise_wts = FALSE#,     # sum (not average) the mean-process over the SA1s
  # known_sigma2fs = 0.2127   # provide estimate of the fine-scale variance
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
# 0.901 ## using SA1s as training data  
# 0.908 ## using mixture of SA1s and SA2s as training data

write.csv(data.frame(coverage = coverage), 
          "./results/Sydney_SA1_coverage.csv", 
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
  S, pred, 
  map_layer = SA2_bg,            # optional layer to put below the plotting geom
  colour = "black", size = 0.025 # optional arguments to plotting geom via ...
  ) 

plots <- lapply(plots, function(gg) gg + lab1 + lab2)

ggsave( 
  ggarrange(plots$p_prob, plots$interval90_prob, nrow = 1, legend = "top"),
  filename = "Sydney_SA1_predictions.png", device = "png", width = 9.5, height = 4.4,
  path = "./img/"
)


# ---- SA3 predictions ----

RNGversion("3.6.0"); set.seed(1)
pred <- predict(S, newdata = SA3_NSW_sub)
plots <- plot(S, pred, colour = "black")
plots <- lapply(plots, function(gg) gg + lab1 + lab2)

ggsave( 
    ggarrange(plots$p_mu, 
              plots$interval90_mu + 
                scale_fill_distiller(palette="BrBG", breaks = c(100, 200, 300)), 
              nrow = 1, legend = "top"),
  filename = "Sydney_SA3_predictions.png", device = "png", width = 9.5, height = 4.4,
  path = "./img/"
)

options(warn = defaultW)

