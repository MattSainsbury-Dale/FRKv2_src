# About:
#   This file conducts the example in Section 3.1 of simulated, spatial, point-
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

suppressMessages({
library("FRK")
library("ggplot2")
library("sp")
library("dplyr")
library("ggpubr")
library("RandomFields")

source("./scripts/Plotting_helpers/Plotting_helpers.R")
source("./scripts/Utility_fns.R")
})

# Define how we wish to simulate the latent process Y
simulation_method <-  "trig_field" # "model" 

## Define a grid of BAUs
x.seq <- y.seq <- seq(0, 1, length.out = 100)
BAUs <- expand.grid(x = x.seq, y = y.seq)
coordinates(BAUs) = ~ x + y
gridded(BAUs) <- TRUE

BAUs_df <- coordinates(BAUs) %>% as.data.frame

## Simulate the true process, Y, over the BAUs
RNGversion("3.6.0"); set.seed(2020)
if (simulation_method == "trig_field") {
  ## Define some complicated smooth trigonometric field
  f <- function(x, y, a = 0, b = 0, l = 1) exp(-l * sqrt((x - a)^2 + (y - b)^2))
  logistic <- function(x, L = 1, k=1, x0 = 0) L / (1 + exp(-k * (x - x0)))
  smooth_Y_process <- function(x, y) {
    a <- 4 + 2 * sin(5 * x) + 2 * cos(4 * y) + 2 * sin(3 * x * y)+
      sin(7 * x) + cos(9 * y) +
      sin(12 * x) + cos(17 * y) + sin(14 * x) * cos(16 * y)  +
      sin(23 * x) + cos(22 * y) + sin(24 * x) * cos(26 * y) +
      x + y + x^2 + y^2 + x * y + 
      sin(10 * pi * x * y) + cos(20 * x * y - 3) + sin(40 * x * y)
    
    logistic(a, x0 = 2, L = 6, k = 0.35)
  }
  
  BAUs_df <- BAUs_df %>%
    mutate(Y = smooth_Y_process(x, y) + rnorm(length(BAUs), mean = 0, sd = 0.2)) 
  
} else if (simulation_method == "model") {
  
  model <- RMexp(var = 0.1, scale = 0.2) + # exponential covariance function
    # RMnugget(var = 0.01) + # nugget
    RMtrend(mean = 4) # and mean
  
  simu <- RFsimulate(model, x = x.seq, y = y.seq)
  BAUs_df$Y <- simu$variable1
}

## Compute the true mean process, mu, over the BAUs
BAUs_df <- mutate(BAUs_df, mu = exp(Y))

## Subsample n of the BAUs to act as observation locations, and simulate data
n <- 750
Poisson_simulated <- sample_n(BAUs_df, n) %>% 
  dplyr::mutate(Z = rpois(n, lambda = mu)) %>%
  dplyr::select(x, y, Z)

# save(Poisson_simulated, file = "~/FRK/data/Poisson_simulated.rda")

## scalar matrix for fine scale variation, and converts to SpatialPixelsDF
BAUs$fs <- rep(1, length(BAUs)) 

## Convert Poisson_simulated to Spatial* obejct
coordinates(Poisson_simulated) <- ~ x + y

# ---- FRK ----

## Predictive performance with changing number of basis functions: Fit the SRE 
## object using 1, 2, and 3 resolutions of basis functions.
max_nres <- 3
pred_list <- S_list <- timings <- list()
for (i in 1:max_nres) {
  timings[[i]] <- system.time({
    S_list[[i]]    <- FRK(f = Z ~ 1, data = list(Poisson_simulated), 
                          nres = i, BAUs = BAUs, 
                          response = "poisson", 
                          link = "log", 
                          # manually set these arguments to reduce console output:
                          K_type = "precision", method = "TMB", est_error = FALSE) 
    RNGversion("3.6.0"); set.seed(1)
    pred_list[[i]] <- predict(S_list[[i]], type = c("link", "mean"))
  })
}


## Predictions, uncertainty, and data
plot_list <- plot(S_list[[max_nres]], pred_list[[max_nres]]$newdata, 
                  labels_from_coordnames = FALSE)
plot_list <- c(plot_spatial_or_ST(Poisson_simulated, "Z", 
                                  labels_from_coordnames = FALSE), 
               plot_list)

## True mean process
plot_list$mu_true <-  ggplot(BAUs_df) + geom_tile(aes(x, y, fill = mu)) + 
  scale_fill_distiller(palette = "Spectral") + 
  labs(fill = bquote(bold("\U03BC"))) +
  labs(x = expression(s[1]), y = expression(s[2])) +
  theme_bw() + coord_fixed()

## Set the title of each plot using the legend labels, 
# and then remove the legend labels. Also adjust the axis breaks, and 
# increase the font size
plot_list <- lapply(
  plot_list,
  function(gg) {
    gg + labs(title = gg$labels$fill) + 
      labs(fill = "", colour = "") + 
      scale_x_continuous(breaks=c(0.25, 0.75), expand = c(0, 0)) + 
      scale_y_continuous(breaks=c(0.25, 0.75), expand = c(0, 0)) + 
      theme(axis.text = element_text(size = 16),
            axis.title = element_text(size = 19), 
            legend.text = element_text(size = 16), 
            plot.title = element_text(hjust = 0.5, size = 19))
  })

## Increase legend width for plots that will have the legend on the side. 
## Also shift the title further to the right, so it looks better centred. 
for (i in c("p_Y", "interval90_Y", "p_mu", "interval90_mu")) {
  plot_list[[i]] <- plot_list[[i]] + 
    theme(legend.key.width = unit(1.1, 'cm'), 
          plot.title = element_text(hjust = 0.66, size = 19))
}

## Increase legend height for plots that will have vertical legends. 
## Also shift the title further to the right, so it looks better centred. 
for (i in c("Z", "mu_true")) {
  plot_list[[i]] <- plot_list[[i]] + 
    theme(legend.key.height = unit(0.8, 'cm'))
}


## Adjust the scale of the mean prediction, so it is the same as the data scale
data_lims  <- range(pred_list[[max_nres]]$newdata$p_mu, Poisson_simulated$Z, BAUs_df$mu)

plot_list$p_mu    <- plot_list$p_mu %>% change_legend_limits(limits = data_lims) 
plot_list$mu_true <- plot_list$mu_true %>% change_legend_limits(limits = data_lims) 
plot_list$Z       <- plot_list$Z %>% change_legend_limits(limits = data_lims, aesthetic = "colour") 




# Adjust the breaks
if (simulation_method == "trig_field") {

  breaks_data <- c(100, 300, 500)
  plot_list$p_mu    <- plot_list$p_mu %>% change_legend_breaks(breaks_data) 
  plot_list$mu_true <- plot_list$mu_true %>% change_legend_breaks(breaks_data) 
  plot_list$Z       <- plot_list$Z %>% change_legend_breaks(breaks_data, aesthetic = "colour") 
  

  plot_list$Z <- plot_list$Z  +
    labs(title = bquote(bold("Z")))
    scale_colour_distiller(palette = "Spectral", name = "", breaks = breaks_data)

    suppressMessages(
  plot_list$interval90_mu <- plot_list$interval90_mu +
    scale_fill_distiller(palette = "BrBG", name = "", breaks = c(100, 250, 400))
  )
}

plot_list$Z <- plot_list$Z +  labs(title = bquote(bold("Z"))) 
figure <- ggarrange(plot_list$mu_true, plot_list$Z, 
                    nrow = 1, common.legend = TRUE, legend = "right", align = "hv")

ggsave(
  figure,
  filename = "3_1_Poisson_sim_true_process_and_data.png", 
  path = "./results", device = "png", width = 10, height = 4
)

## Remove y-axis labels/ticks for all but the left-most panel
interior_plot <- function(gg) {
  gg + rremove("ylab") + rremove("xlab") + rremove("y.text") + rremove("y.ticks")
}

exterior_plot <- function(gg) {
  gg + rremove("ylab") + rremove("xlab")
}


figure <- ggarrange(plot_list$p_Y %>% exterior_plot,
                    plot_list$interval90_Y %>% interior_plot, 
                    plot_list$p_mu %>% interior_plot, 
                    plot_list$interval90_mu %>% interior_plot, 
                    nrow = 1, legend = "top", align = "hv") %>% 
  annotate_figure(left = text_grob(bquote(s[2]), rot = 90, vjust = 1, hjust = 2, size = 20),
                  bottom = text_grob(bquote(s[1]), size = 20, vjust = -1, hjust = -0.5))


ggsave(
  figure,
  filename = "3_1_Poisson_sim.png", 
  path = "./results", device = "png", width = 14, height = 5.2
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
unobsidx <- unobserved_BAUs(S_list[[1]]) # doesn't matter which SRE object we use 

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

write.csv(diagnostics, file = "./results/3_3_Poisson_nres_comparison.csv")

save_html_table(
  diagnostics,
  file = "results/3_3_Poisson_nres_comparison.html", 
  caption = "Poisson comparison using multiple resolutions of basis functions"
)