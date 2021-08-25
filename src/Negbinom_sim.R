## Toy example for spatial COS.
library("FRK")
library("sp")
library("ggplot2")
library("ggpubr")
library("dplyr")
library("raster")

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
xbreaks <- scale_x_continuous(breaks=c(0.25, 0.75), expand = c(0, 0))
ybreaks <- scale_y_continuous(breaks=c(0.25, 0.75), expand = c(0, 0)) 

## Legend breaks
breaks_prob = c(0.25, 0.5, 0.75)
breaks_mu_BAU = c(50, 125, 200)
breaks_data = c(2500, 6000, 9500)

## Now create the mean process over the BAUs. 
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

g_prob_BAU <-  plot_spatial_or_ST(BAUs, "prob", labels_from_coordnames = F)[[1]] + 
  scale_fill_distiller(palette = "Spectral", breaks = breaks_prob) + 
  labs(title = bquote(bold("\U03C0")), fill = "") + 
  xbreaks + ybreaks 

## We need the size parameter at the BAU level
BAUs$k_BAU <- rep(50, length(BAUs))

## Mean over the BAUs
BAUs$mu <- BAUs$k_BAU * (1 / BAUs$prob - 1)

g_mu_BAU <- plot_spatial_or_ST(BAUs, "mu", labels_from_coordnames = F)[[1]] + 
  scale_fill_distiller(palette = "Spectral", breaks = breaks_mu_BAU) + 
  labs(title = bquote(bold("\U03BC")), fill = "") + 
  xbreaks + ybreaks 

## Now aggregate the mean process over the data supports
BAUs_as_points <- SpatialPointsDataFrame(coordinates(BAUs), 
                                         BAUs@data, 
                                         proj4string = CRS(proj4string(BAUs)))
obs <- SpatialPolygonsDataFrame(obs, 
                                over(obs, BAUs_as_points, fn = sum))

## Use a subsample of the data supports for model fitting. 
obsidx <- sample(1:length(obs), round(length(obs) * 0.8), replace = FALSE)
obs <- obs[obsidx, ]

## Plot the aggregated mean process
g_mu <- plot_spatial_or_ST(obs, "mu", colour = "black", size = 0.1, 
                           labels_from_coordnames = F)[[1]] + 
  scale_fill_distiller(palette = "Spectral", breaks = breaks_data) + 
  labs(title = bquote(bold("\U03BC")[Z]), fill = "") + 
  xbreaks + ybreaks
  

## Now simulate data over the data supports. 
obs$Z <- rnbinom(n = length(obs), size = obs$k_BAU, mu = obs$mu)


zdf <- obs

## Training data
g_Z_training <- plot_spatial_or_ST(zdf, "Z", colour = "black", size = 0.1, 
                                   labels_from_coordnames = F)[[1]] + 
  xbreaks + ybreaks + 
  scale_fill_distiller(palette = "Spectral", breaks = breaks_data) + 
  labs(title = bquote(bold("Z")), fill = "")


set_title_from_fill_legend <- function(gg) {
  gg + labs(title = gg$labels$fill) + labs(fill = "")
}

change_font_size <- function(gg, size = 11) {
  gg + theme(axis.text = element_text(size = 16),
             axis.title = element_text(size = 19), 
             legend.text = element_text(size = 16), 
             plot.title = element_text(hjust = 0.5, size = 19))
}

change_legend_width <- function(gg, width = 1.1) {
  gg + theme(legend.key.width = unit(width, 'cm'))
}



data_plots <- list(g_prob_BAU, g_mu_BAU, g_mu, g_Z_training)
data_plots <- lapply(data_plots, change_font_size)
data_plots <- lapply(data_plots, change_legend_width)

## Remove y-axis labels/ticks for all but the left-most panel
interior_plot <- function(gg) {
  gg + rremove("ylab") + rremove("xlab") + rremove("y.text") + rremove("y.ticks")
}
exterior_plot <- function(gg) {
  gg + rremove("ylab") + rremove("xlab")
}

create_figure_one_row_of_plots <- function(plot_list) {
  plot_list[[1]] = plot_list[[1]] %>% exterior_plot
  for (i in 2:length(plot_list)) {
    plot_list[[i]] = plot_list[[i]] %>% interior_plot
  }
  
  figure <- ggarrange(plotlist = plot_list, 
                      nrow = 1, legend = "top", align = "hv") %>% 
    annotate_figure(left   = text_grob(bquote(s[2]), size = 20, vjust = 1,  hjust = 2, rot = 90),
                    bottom = text_grob(bquote(s[1]), size = 20, vjust = -1, hjust = -0.5))
} 


figure <- create_figure_one_row_of_plots(data_plots)

ggsave(figure,
       filename = "Negbinom_sim_data.png", device = "png", 
       width = 14, height = 5.1,
       path = "./img"
)



# ---- Prep BAUs and zdf for FRK ----

BAUs$fs <- 1 # scalar variance matrix for fine-scale random effects

## FRK expects the size parameter at the data level to be labelled 'k_Z'
names(zdf)[names(zdf) == "k_BAU"] <- "k_Z"


## Remove the true processes from zdf and BAUs, and Z from BAUs
## (so we could not possibly use them in fitting). 
## Save the true mean for validation later. 
mu_true <- BAUs$mu
BAUs@data[, c("x", "y", "Y", "prob", "mu", "Z")] <- NULL
zdf@data[, c("x", "y", "Y", "prob", "mu")] <- NULL
## (keep BAUs$mu for model validation)

# ---- FRK ----

S <- FRK(f = Z ~ 1, data = list(zdf), BAUs = BAUs, 
         response = response, link = link, 
         normalise_wts = FALSE, method = "TMB")

pred <- predict(S)


# ---- Plotting results ----

plot_list <- plot(S, pred$newdata, labels_from_coordnames = F)

plot_list <- lapply(plot_list, function(gg) gg + xbreaks + ybreaks)
plot_list <- lapply(plot_list, change_font_size)
plot_list <- lapply(plot_list, change_legend_width)
plot_list <- lapply(plot_list, set_title_from_fill_legend)

plot_list$p_prob <- plot_list$p_prob + scale_fill_distiller(palette = "Spectral", breaks = breaks_prob, name = "") 
plot_list$interval90_prob <- plot_list$interval90_prob + scale_fill_distiller(palette = "BrBG", breaks = c(0.1, 0.15, 0.2), name = "")
plot_list$p_mu <- plot_list$p_mu + scale_fill_distiller(palette = "Spectral", breaks = breaks_mu_BAU, name = "")
plot_list$interval90_mu <- plot_list$interval90_mu + scale_fill_distiller(palette = "BrBG", breaks = c(50, 100, 150), name = "")

figure <- create_figure_one_row_of_plots(list(plot_list$p_prob, 
                                              plot_list$interval90_prob, 
                                              plot_list$p_mu, 
                                              plot_list$interval90_mu))

ggsave( 
  figure,
  filename = "Negbinom_sim_BAU_predictions.png", device = "png", 
  width = 14, height = 5.2,
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
    true_in_pred_interval90 = mu_percentile_5 <= true & true <= mu_percentile_95,
    absolute_percentage_error = abs(p_mu - true)/true
  ) %>%
  summarise(
    coverage_90 = mean(true_in_pred_interval90),
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

## NB: Cannot make predictions of pi(.) over arbitrary polygons, only the mean process
pred_over_polygons <- predict(S, newdata = newdata)

plot_list <- plot(S, pred_over_polygons$newdata, colour = "black", labels_from_coordnames = F)

## Hack to make the correct x-axis and y-axis breaks appear
invisible_point <- geom_point(data = data.frame(x = c(0, 1), y = c(1, 1)), alpha = 0)
plot_list <- lapply(plot_list, function(gg) gg + invisible_point)

plot_list <- lapply(plot_list, function(gg) gg + xbreaks + ybreaks)
plot_list <- lapply(plot_list, change_font_size)
plot_list <- lapply(plot_list, set_title_from_fill_legend)
plot_list <- lapply(plot_list, function(gg) gg + theme(legend.key.height = unit(1.1, "cm")))

plot_list$p_mu <- plot_list$p_mu + 
  scale_fill_distiller(palette = "Spectral", breaks = c(8000, 12000, 16000))

ggarrange(plot_list$p_mu, plot_list$interval90_mu, nrow = 1, align = "hv")
ggarrange(plot_list$p_prob, plot_list$interval90_prob, nrow = 1, align = "hv")

ggsave( 
  ggarrange(plot_list$p_mu, plot_list$interval90_mu, nrow = 1, align = "hv"),
  filename = "Negbinom_sim_arbitrary_polygon_predictions.png", device = "png", 
  width = 11, height = 4.8,
  path = "./img"
)


