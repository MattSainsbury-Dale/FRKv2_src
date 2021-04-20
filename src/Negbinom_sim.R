## Toy example for spatial COS.
library("FRK")
library("sp")
library("ggplot2")
library("ggpubr")
library("dplyr")

## set the seed once at the start of the script. Results should be reproducible.
RNGversion("3.6.0")
set.seed(1)

## Some helper functions to construct SpatialPolygons used as data supports, 
## BAUs and arbitrary prediction polygons.
source("./src/SpatialPolygon_fns.R")


# ---- Construct BAUs, data supports, and arbitrary regions to predict over ----

BAUs <- makeGrid(n = 64)
obs <- makeGrid(n = 8)

## Randomly generate some non-regular convex polygons
## NB: use same area so that comparison between polygons is fair
coords <- mapply(convex.poly, nSides = 3:6, area = 0.07) 

## Shift each polygon by "a"
coords <- mapply(shift.coordinates, coords, 
                 a = list(c(0.25, 0.28), 
                          c(0.4, 0.7), 
                          c(0.8, 0.7), 
                          c(0.8, 0.2)))

newdata <- coords_to_SpatialPolygons(coords)
coordnames(BAUs) <- coordnames(obs) <- coordnames(newdata) <- c("x", "y")

# ---- Use BAUs for some data observations ----

## Select some data supports to be replaced
rmidx <- sample(1:length(obs), ceiling(length(obs)/10), replace = FALSE)

## Convert BAUs to points, otherwise BAUs in neighbouring data supports are selected too, 
## because they share lines, which is classified as overlapping by over()
BAUs_as_points <- SpatialPoints(BAUs)
## Select BAUs corresponding to the selected data supports
observed_BAU_idx <- which(!is.na(over(BAUs_as_points, obs[rmidx])))
## Combine selected BAUs and remaining data supports
obs <- raster::bind(BAUs[observed_BAU_idx], obs[-rmidx])


# ---- Define true processes Y(.), pi(.), and mu(.), and simulate data ----

response <- "negative-binomial" # "Gaussian"
link <- "logit" # "identity"

## Define the smooth process (responsible for X*beta + S*eta)
smooth_Y_process <- function(x, y) sin(11 * x) + cos( 10 * y)

## Simulate true process Y
sigma2fs <- 0.05
BAUs <- SpatialPixels(SpatialPoints(BAUs))
## We need to convert the BAUs to SpatialPixelsDataFrame: just do this by
## setting fs = 1 (indicating the fine-scale variation is homoscedastic)
BAUs$fs <- rep(1, length(BAUs)) 
BAUs@data <- as.data.frame(coordinates(BAUs))
BAUs$Y <- smooth_Y_process(BAUs$x, BAUs$y) + rnorm(length(BAUs), sd = sqrt(sigma2fs))

## Common axis breaks
xbreaks <- scale_x_continuous(breaks=c(0, 0.5, 1)) 
ybreaks <- scale_y_continuous(breaks=c(0, 0.5, 1))

## Now create the mean process over the BAUs. 
if (link == "logit") {
  
  ## First construct the probability process
  BAUs$prob <- plogis(BAUs$Y)
  
  ## For visually appealing data, we don't want the probability process 
  ## to be extremely close to 0 or 1, so lets squash it to be between, say, 
  ## 0.1 and 0.9. We'll do this using the general logistic function shifted 
  ## by a constant. As the input, p_O, is already restricted to be in [0, 1],
  ## we will set the curves centre to be 0.5.
  general_logistic <- function(x, L = 1, k = 1, x0 = 0, a = 0){
    L / (1 + exp(-k * (x - x0))) + a
  }
  
  BAUs$prob <- general_logistic(BAUs$prob, L = 0.95, a = 0.05, x0 = 0.5, k = 4)
  
  g_prob_BAU <-  plot_spatial_or_ST(BAUs, "prob")[[1]] + 
    labs(fill = expression(pi("\U00B7"))) + 
    xbreaks + ybreaks 
  
  ## We need the size parameter at the BAU level
  BAUs$k_BAU <- rep(50, length(BAUs))
  
  if (response == "negative-binomial") 
    BAUs$mu <- BAUs$k_BAU * (1 / BAUs$prob - 1)
  
} else if (link == "identity") {
  BAUs$mu <- BAUs$Y
}

g_mu_BAU <- plot_spatial_or_ST(BAUs, "mu")[[1]] + 
  labs(fill = expression(mu("\U00B7"))) + 
  xbreaks + ybreaks 

## Now aggregate the mean process over the data supports
BAUs_as_points <- SpatialPointsDataFrame(coordinates(BAUs), 
                                         BAUs@data, 
                                         proj4string = CRS(proj4string(BAUs)))
obs <- SpatialPolygonsDataFrame(obs, 
                                over(obs, BAUs_as_points, fn = sum))

## breaks for the data-support scale
data_breaks <- c(3000, 6000, 9000)

## Plot the aggregated mean process
g_mu <- plot_spatial_or_ST(obs, "mu", colour = "black", size = 0.1)[[1]] + 
  labs(fill = expression(mu("\U00B7"))) + 
  xbreaks + ybreaks + 
  scale_fill_distiller(palette = "Spectral", breaks = data_breaks)

## Now simulate data over the data supports. 
if (response == "negative-binomial") {
  obs$Z <- rnbinom(n = length(obs), size = obs$k_BAU, mu = obs$mu)
} else if (response == "Gaussian") {
  sigma2e <- 3^2
  obs$Z <- rnorm(n = length(obs), mean = obs$mu, sd = sqrt(sigma2e))
}

## Use a subsample of the data for model fitting
obsidx <- sample(1:length(obs), round(length(obs) * 0.8), replace = FALSE)
zdf <- obs[obsidx, ]

## Training data
g_Z_training <- plot_spatial_or_ST(zdf, "Z", colour = "black", size = 0.1)[[1]] + 
  xbreaks + ybreaks + scale_fill_distiller(palette = "Spectral", breaks = data_breaks)

if (response == "negative-binomial") {
  ggsave( 
    ggarrange(g_prob_BAU, g_mu_BAU, g_mu, g_Z_training, 
              nrow = 1, align = "hv", legend = "top"),
    filename = "Negbinom_sim_data.png", device = "png", 
    width = 13, height = 4,
    path = "./img"
  )
} else if (response == "Gaussian") {
  ggarrange(g_mu_BAU, g_mu, g_Z_training, nrow = 1, align = "hv", legend = "top") 
}


# ---- Prep BAUs and zdf for FRK ----

BAUs$fs <- 1 # scalar variance matrix for fine-scale random effects

## FRK expects the size parameter at the data level to be labelled 'k_Z'
if (response == "negative-binomial") {
  names(zdf)[names(zdf) == "k_BAU"] <- "k_Z"
}

## Remove the true processes from zdf and BAUs, and Z from BAUs
## (so we could not possibly use them in fitting). 
## Save the true mean for validation later. 
mu_true <- BAUs$mu
BAUs@data[, c("x", "y", "Y", "prob", "mu", "Z")] <- NULL
zdf@data[, c("x", "y", "Y", "prob", "mu")] <- NULL
## (keep BAUs$mu for model validation)

# ---- FRK ----

est_error <- FALSE
if (response == "Gaussian" && !est_error) zdf$std <- 1

S <- FRK(f = Z ~ 1, data = list(zdf), BAUs = BAUs, 
         response = response, link = link, 
         normalise_wts = FALSE, method = "TMB", 
         est_error = est_error)

pred <- predict(S)


# ---- Plotting results ----

plot_list <- plot(S, pred)

plot_list <- lapply(
  plot_list,
  function(gg) gg + xbreaks + ybreaks)

ggsave( 
  ggarrange(plot_list$p_prob, plot_list$interval_90_prob, 
            plot_list$p_mu, plot_list$interval_90_mu, 
            nrow = 1, legend = "top"),
  filename = "Negbinom_sim_BAU_predictions.png", device = "png", 
  width = 13, height = 4,
  path = "./img"
)



# ---- Assess prediction and UQ at BAU level ----

## Compare BAUs$mu with the predictions (i.e., compute RMSPE), and the coverage
## via the percentiles.
pred$newdata$true <- mu_true

## RMSPE
RMSPE <- function(z,pred) sqrt(mean((z - pred)^2))
RMSPE(mu_true, pred$newdata$p_mu)  

## Coverage and MAPE
diagnostics <- pred$newdata@data %>% 
  mutate(
    true_in_pred_interval_90 = mu_percentile_5 <= true & true <= mu_percentile_95,
    absolute_percentage_error = abs(p_mu - true)/true
  ) %>%
  summarise(
    coverage_90 = mean(true_in_pred_interval_90),
    MAPE = mean(absolute_percentage_error) * 100, 
    RMSPE = RMSPE(true, p_mu)  
  ) 

## CRPS
diagnostics <- cbind(
  diagnostics,
  CRPS = mean(scoringRules::crps_sample(y = mu_true, dat = pred$MC$mu_samples))
)

write.csv(diagnostics, file = "./results/Negbinom_sim.csv", row.names = FALSE)


# ---- Predict over polygons ----

##NB: Cannot make predictions of pi(.) over arbitrary polygons, only the mean process
pred_over_polygons <- predict(S, newdata = newdata)

plot_list <- plot(S, pred_over_polygons, colour = "black")

ggsave( 
  ggarrange(plot_list$p_mu, plot_list$interval_90_mu, nrow = 1, align = "hv"),
  filename = "Negbinom_sim_arbitrary_polygon_predictions.png", device = "png", 
  width = 10, height = 3,
  path = "./img"
)
