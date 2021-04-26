# File: resolution_comparison.R
# Author: Matthew Sainsbury-Dale
# About:
#   This file conducts the example in Section 3.1 of simulated spatial, point-
#   referenced, Poisson-distributed data, including producing plots of the true
#   mean process mu(.), the data Z, and predictions and associated uncertainties 
#   of the latent process Y(.) and mean process mu(.). 
#   It also compares the performance of FRK when using basis functions of 1, 2, 
#   and 3 resolutions. It does so using prediction and uncertainty plots at each 
#   resolution, as well as computing the diagnostics based on the true process 
#   values. The way in which we generate the true mean process mu(.) is very 
#   important, and the steps are as follows;
#     1. Generate a fine grid to act as BAUs
#     2. Generate a true process (with fine-scale variation) over those BAUs, 
#        defined using trignometric and exponential functions
#     3. Subsample some BAUs to act as observations (call them Z)
#     4. Run FRK using Z, then compute diagnostics using predictions and 
#        previously computed true process value.

library("FRK")
library("ggplot2")
library("sp")
library("dplyr")
library("ggpubr")

## Define some complicated smooth process 
f <- function(x, y, a = 0, b = 0, l = 1) exp(-l * sqrt((x - a)^2 + (y - b)^2))
logistic <- function(x, L = 1, k=1, x0 = 0) L / (1 + exp(-k * (x - x0)))
smooth_Y_process <- function(x, y) {
  a <- 2 * sin(5 * x) + 2 * cos(4 * y) + 2 * sin(3 * x * y)+
    sin(7 * x) + cos(9 * y) +
    sin(12 * x) + cos(17 * y) + sin(14 * x) * cos(16 * y)  +
    sin(23 * x) + cos(22 * y) + sin(24 * x) * cos(26 * y) +
    x + y + x^2 + y^2 + x * y + 4 +
    sin(10 * pi * x * y) + cos(20 * x * y - 3) + sin(40 * x * y)
  
  logistic(a, x0 = 2, L = 6, k = 0.35)
}

## Define a grid of BAUs
BAUs <- expand.grid(x = seq(0, 1, length.out = 100),
                    y = seq(0, 1, length.out = 100))
coordinates(BAUs) = ~ x + y
gridded(BAUs) <- TRUE

## Generate the true process Y(.) (with fine-scale variation) over the BAUs,
## and the true mean process, mu(.)
BAUs_df <- coordinates(BAUs) %>%
  as.data.frame() %>%
  mutate(Y = smooth_Y_process(x, y) + rnorm(length(BAUs), mean = 0, sd = 0.2)) %>%
  mutate(mu = exp(Y))

## Subsample n of the BAUs to act as observation locations, and simulate data
n <- 750
Poisson_simulated <- sample_n(BAUs_df, n) %>% 
  mutate(Z = rpois(n, lambda = mu)) %>%
  select(x, y, Z)

# save(Poisson_simulated, file = "~/FRK/data/Poisson_simulated.rda")

## scalar matrix for fine scale variation, and converts to SpatialPixelsDF
BAUs$fs <- rep(1, length(BAUs)) 

## Convert Poisson_simulated to Spatial* obejct
coordinates(Poisson_simulated) <- ~ x + y



# ---- Predictive performance with changing number of basis functions ----

max_nres <- 3
pred_list <- S_list <- timings <- list()
for (i in 1:max_nres) {
  S_list[[i]]    <- FRK(f = Z ~ 1, data = list(Poisson_simulated), 
              nres = i, BAUs = BAUs, 
              response = "poisson", link = "log")
  RNGversion("3.6.0"); set.seed(1)
  pred_list[[i]] <- predict(S_list[[i]], type = c("link", "mean"))
}

## Predictions, uncertainty, and data
plot_list <- plot(S_list[[3]], pred_list[[3]], Poisson_simulated)

## True mean process
plot_list$mu_true <-  ggplot(BAUs_df) + geom_tile(aes(x, y, fill = mu)) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill = bquote("\U03BC(\U00B7)")) +
  theme_bw() + coord_fixed()

## Adjust the breaks and remove the fill name
plot_list <- lapply(
  plot_list, 
  function(gg) gg + scale_x_continuous(breaks=c(0, 0.5, 1)) + 
    scale_y_continuous(breaks=c(0, 0.5, 1)))

## Adjust the scale of the mean prediction, so it is the same as the data scale
data_scale <- range(pred_list[[3]]$newdata$p_mu, Poisson_simulated$Z)
plot_list$p_mu <- plot_list$p_mu +  scale_fill_distiller(palette = "Spectral", 
                                                         limits = data_scale, 
                                                         breaks = c(100, 300, 500))



ggsave(
  ggarrange(plot_list$mu_true, plot_list$Z, nrow = 1),
  filename = "Poisson_sim_true_process_and_data.png", 
  path = "./img", device = "png", width = 10, height = 4
)

ggsave(
  ggarrange(plot_list$p_Y, plot_list$interval90_Y, 
            plot_list$p_mu, plot_list$interval90_mu, 
            nrow = 1, legend = "top"),
  filename = "Poisson_sim.png", 
  path = "./img", device = "png", width = 14, height = 4
)


## Diagnostic functions
.RMSPE <- function(pred, true) sqrt(mean((pred - true)^2))
.Coverage <- function(lower, upper, true) ((true > lower) & (true < upper))
.IS90 <- function(lower, upper, true, alpha = 0.1) {
  (upper - lower) + 2/alpha * (lower - true) * (true < lower) +
    2/alpha * (true - upper) * (true > upper)
}
.diagnostic_stats <- function(lower, upper, pred, samples, true) {
  data.frame(RMSPE = .RMSPE(true, pred),
             CRPS = mean(scoringRules::crps_sample(y = true, dat = samples)),
             IS90 = mean(.IS90(lower, upper, true)),
             Coverage = mean(.Coverage(lower, upper, true)))
}

## Only consider out-of-sample locations
unobsidx <- unobserved_BAUs(S) 

## Compute the diagnostics
diagnostics <- lapply(pred_list, function(pred) {
  .diagnostic_stats(
    lower = pred$newdata$mu_percentile_5[unobsidx],
    upper = pred$newdata$mu_percentile_95[unobsidx],
    pred = pred$newdata$p_mu[unobsidx],
    samples = pred$MC$mu_samples[unobsidx, ],
    true = BAUs_df$mu[unobsidx])
}) 

diagnostics <- do.call(rbind.data.frame, diagnostics)

diagnostics <- cbind(diagnostics, 
                     run_time = sapply(timings, function(x) x["elapsed"]))

write.csv(diagnostics, file = "./results/Poisson_nres_comparison.csv")

# ## Extract prediction and UQ from each resolution, and convert to single dataframe.
# newdata_list <- lapply(pred_list, function(x) x$newdata)
# newdata <- do.call("rbind", newdata_list)
# newdata$nres <- rep(paste("nres =", 1:nrow(diagnostics)), 
#                     each = length(newdata_list[[1]]))
# 
# newdata@data <- newdata@data %>% cbind(coordinates(newdata))
# 
# g_basic <- ggplot(newdata@data, aes(x, y)) + 
#   theme_bw() + coord_fixed() + 
#   scale_x_continuous(breaks=c(0, 0.5, 1)) + 
#   scale_y_continuous(breaks=c(0, 0.5, 1))
# 
# g_mu_pred <- g_basic + geom_tile(aes(fill = p_mu)) +
#   scale_fill_distiller(palette = "Spectral") +
#   labs(fill = expression(widehat(p)[mu]["|"][bold(Z)])) + 
#   facet_wrap(~nres, nrow = 1)
# 
# g_mu_UQ <- g_basic + geom_tile(aes(fill = mu_percentile_95 - mu_percentile_5)) +
#   scale_fill_distiller(palette = "BrBG") + 
#   labs(fill = "90% predictive\ninterval width\nof the mean process") + 
#   facet_wrap(~nres, nrow = 1)
# 
# ggarrange(g_mu_pred, g_mu_UQ, nrow = 2, align = "hv") # Export: 13.1 x 7.9 in
