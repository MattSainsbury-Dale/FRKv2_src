suppressMessages({
library("FRK")
library("sp")
library("ggplot2")
library("ggpubr")
library("dplyr")
library("raster")
# library("RandomFields")
})

source("Code/Utility_fns.R")

## Use very-low-dimensional representations of the models to establish that the code works? 
quick <- check_quick()
nres  <- if (quick) 2 else 3

# Define how we wish to simulate the latent Y process
simulation_method <- "trig_field"  # model

RNGversion("3.6.0"); set.seed(1)


# ---- Helper functions ----

## Some helper functions to construct SpatialPolygons used as data supports, 
## BAUs and arbitrary prediction polygons, and plotting helper functions.

## Lower left-hand vertex is (a, b), side length s
makeSquare <- function(a, b, s) {
  
  ## Start with lower left-hand vertex then continue clockwise
  vertex_corners <- rbind(c(a, b), 
                          c(a, b + s), 
                          c(a + s, b + s), 
                          c(a + s, b))
  
  square <- Polygon(vertex_corners)
  
  return(square)
}

## s is the side length of the grid squares, 
## n is the number of grid squares in each dimension,
makeGrid <- function (n, range0 = 1) { 
  
  ## range0 is the desired extent of the grid.
  ## Define s, the side length of the squares, so that s*n = range = range0:
  s <- range0 / n
  
  squares <- list()
  i <- 1
  for (b in seq(0, (n-1)*s, by = s)) {
    for (a in seq(0, (n-1)*s, by = s)) {
      square <- Polygons(list(makeSquare(a = a, b = b, s = s)), 
                         paste0("square", i))
      squares[[i]] <- square
      i = i + 1
    }
  }
  
  square_grid <- SpatialPolygons(squares)
  
  return(square_grid)
}


## Shift coordinates by a: a is scalar or length-two vector
shift.coordinates <- function(coords, a) {
  
  if (length(a) == 1) {
    coords <- coords + a
  } else if (length(a) == 2) {
    coords <- coords + rep(a, each = nrow(coords))
  } else 
    stop("a should be a scalar or a vector of length 2.")
}


## Credit to this function goes to:
# https://stackoverflow.com/a/4860337
convex.poly <- function(nSides, area) 
{
  # Find the radius of the circumscribed circle, and the angle of each point if this was a regular polygon
  radius <- sqrt((2*area)/(nSides*sin((2*pi)/nSides)))
  angle <- (2*pi)/nSides
  
  # Randomize the radii/angles
  radii <- rnorm(nSides, radius, radius/10)
  angles <- rnorm(nSides, angle, angle/10) * 1:nSides
  angles <- sort(angles)
  
  points <- list(x=NULL, y=NULL)
  points$x <- cos(angles) * radii
  points$y <- sin(angles) * radii
  
  # Find the area of the polygon
  m <- matrix(unlist(points), ncol=2)
  m <- rbind(m, m[1,])
  current.area <- 0.5 * (sum(m[1:nSides,1]*m[2:(nSides+1),2]) - sum(m[1:nSides,2]*m[2:(nSides+1),1]))
  
  points$x <- points$x * sqrt(area/current.area)
  points$y <- points$y * sqrt(area/current.area)
  
  ## Append the first set of points to the end, to "complete" the polygon
  points$x <- append(points$x, points$x[1])
  points$y <- append(points$y, points$y[1])
  
  ## Convert to matrix
  points <- do.call(cbind, points)
  
  return (points)
}

## coords: list of matrices containing the coordinates of the polygons
## This function is based on the simple example from vignette("sp")
coords_to_SpatialPolygons <- function(coords) {
  
  ## Convert to Polygon
  Poly_list <- lapply(coords, Polygon)
  
  ## Give each Polygon a name
  names(Poly_list) <- paste0("s", 1:length(Poly_list))
  
  ## Convert to Polygons
  Polygons_list <- lapply(seq_along(Poly_list), function(i) Polygons(list(Poly_list[[i]]), names(Poly_list)[i]))
  
  ## Convert to SpatialPolygons
  newdata <- SpatialPolygons(Polygons_list, 1:length(Polygons_list))
  
  return(newdata)
}



change_font_size <- function(gg) {
  gg + theme(axis.text = element_text(size = 16),
             axis.title = element_text(size = 19), 
             legend.text = element_text(size = 16), 
             plot.title = element_text(hjust = 0.5, size = 19))
}

change_legend_width <- function(gg, width = 1.1) {
  gg + theme(legend.key.width = unit(width, 'cm'))
}

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

## Select some coarse-scale data supports to be replaced with fine-scale data supports
rmidx <- sample(1:length(obs), ceiling(length(obs)/10), replace = FALSE)

## Convert BAUs to points, otherwise BAUs in neighbouring data supports are selected too, 
## because they share lines, which is classified as overlapping by over()
BAUs_as_points <- SpatialPoints(BAUs)
## Select BAUs corresponding to the selected data supports
observed_BAU_idx <- which(!is.na(over(BAUs_as_points, obs[rmidx])))
## Combine selected BAUs and remaining data supports
obs <- raster::bind(BAUs[observed_BAU_idx], obs[-rmidx])


# ---- Define true processes Y(.), pi(.), and mu(.), and simulate data ----

## A list of data plots which we will add to during the simulation procedure.
data_plots <- list()


if (simulation_method == "trig_field") {
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
  
  
} else if (simulation_method == "model") {
  
  RNGversion("3.6.0"); set.seed(1996)
  model <- RMexp(var = 0.6, scale = 1) + # exponential covariance function
    # RMnugget(var = 0.005) + # nugget
    RMtrend(mean = 0) # and mean
  
  ## Simulate true process Y
  BAUs <- SpatialPixels(SpatialPoints(BAUs))
  ## We need to convert the BAUs to SpatialPixelsDataFrame: just do this by
  ## setting fs = 1 (indicating the fine-scale variation is homoscedastic)
  BAUs$fs <- rep(1, length(BAUs)) 
  BAUs@data <- as.data.frame(coordinates(BAUs))
  
  simu <- RFsimulate(model, x = sort(unique(BAUs$x)), y = sort(unique(BAUs$y)))
  BAUs$Y <- simu$variable1
 
  ## Construct the probability process
  BAUs$prob <- plogis(BAUs$Y)
}

## Plot the probability process
data_plots$prob_BAU <-
  plot_spatial_or_ST(BAUs, "prob", labels_from_coordnames = F)[[1]] +  
  labs(title = bquote(bold("\U03C0")), fill = "")

## Mean over the BAUs (requires size parameter at the BAU level)
BAUs$k_BAU <- rep(50, length(BAUs))
BAUs$mu    <- BAUs$k_BAU * (1 / BAUs$prob - 1)

data_plots$mu_BAU <- 
  plot_spatial_or_ST(BAUs, "mu", labels_from_coordnames = F)[[1]] +
  labs(title = bquote(bold("\U03BC")), fill = "")




## Now aggregate the mean process over the data supports
BAUs_as_points <- SpatialPointsDataFrame(coordinates(BAUs), 
                                         BAUs@data, 
                                         proj4string = CRS(proj4string(BAUs)))
obs <- SpatialPolygonsDataFrame(obs, over(obs, BAUs_as_points, fn = sum))

## Use a subsample of the data supports for model fitting
obsidx <- sample(1:length(obs), round(length(obs) * 0.8), replace = FALSE)
obs <- obs[obsidx, ]

## Plot the aggregated mean process
data_plots$mu <- 
  plot_spatial_or_ST(
    obs, "mu", colour = "black", size = 0.1, labels_from_coordnames = F
    )[[1]] +
  labs(title = bquote(bold("\U03BC")[Z]), fill = "") 
  
## Simulate data over the data supports
obs$Z <- rnbinom(n = length(obs), size = obs$k_BAU, mu = obs$mu)

## Training data
data_plots$Z_training <- 
  plot_spatial_or_ST(
    obs, "Z", colour = "black", size = 0.1, labels_from_coordnames = F
    )[[1]] +
  labs(title = bquote(bold("Z")), fill = "")


## Create the data figure. First, define some common plot elements:
xbreaks <- scale_x_continuous(breaks=c(0.25, 0.75), expand = c(0, 0))
ybreaks <- scale_y_continuous(breaks=c(0.25, 0.75), expand = c(0, 0))

                                  
# Adjust the breaks and other plot visuals
if (simulation_method == "trig_field") {
  
  breaks_prob   = c(0.25, 0.5, 0.75)
  breaks_mu_BAU = c(50, 125, 200)
  breaks_data   = c(2500, 6000, 9500)
  breaks <- list(
    prob_BAU = breaks_prob,
    mu_BAU = breaks_mu_BAU,
    mu = breaks_data,
    Z_training = breaks_data
  )
  
  
  data_plots <- lapply(
    names(data_plots), function(i) {
      gg <- data_plots[[i]] 
      gg <- gg %>% change_legend_breaks(breaks[[i]]) 
      gg <- gg + xbreaks + ybreaks
      gg <- gg %>% change_font_size
      gg <- gg %>% change_legend_width
      return(gg)
    }
  )
  
  
} else if (simulation_method == "model") {
  
  data_scale_lims <- range(c(obs$Z, BAUs$mu))
  data_scale_breaks <- c(750, 2000, 3250)
  data_plots$mu <- change_legend_limits(data_plots$mu,  data_scale_lims)
  data_plots$Z_training <- change_legend_limits(data_plots$Z_training,  data_scale_lims)
  data_plots$mu <- change_legend_breaks(data_plots$mu,  data_scale_breaks)
  data_plots$Z_training <- change_legend_breaks(data_plots$Z_training,  data_scale_breaks)
  
  prob_scale_lims <- range(BAUs$prob)
  prob_scale_breaks <- c(0.4, 0.6, 0.8)
  data_plots$prob_BAU <- change_legend_limits(data_plots$prob_BAU,  prob_scale_lims)
  data_plots$prob_BAU <- change_legend_breaks(data_plots$prob_BAU,  prob_scale_breaks)
  
  
  data_plots <- lapply(
    names(data_plots), function(i) {
      gg <- data_plots[[i]] 
      gg <- gg %>% change_legend_breaks(breaks[[i]]) 
      gg <- gg + xbreaks + ybreaks
      gg <- gg %>% change_font_size
      gg <- gg %>% change_legend_width
      return(gg)
    }
  )
}


figure <- create_figure_one_row_of_plots(data_plots)


ggsave(figure,
       filename = "3_2_Negbinom_sim_data.png", device = "png", 
       width = 14, height = 5.1,
       path = "results")


# ---- Prep BAUs and obs for FRK ----

BAUs$fs <- 1 # scalar variance matrix for fine-scale random effects

## FRK expects the size parameter at the data level to be labelled 'k_Z'
names(obs)[names(obs) == "k_BAU"] <- "k_Z"

## Remove the true processes from obs and BAUs, and Z from BAUs so that we 
## could not possibly use them in model fitting. 
## Save the true mean for validation later. 
mu_true <- BAUs$mu
BAUs@data[, c("x", "y", "Y", "prob", "mu", "Z")] <- NULL
obs@data[, c("x", "y", "Y", "prob", "mu")] <- NULL

# ---- Model fitting ----

S <- FRK(f = Z ~ 1, data = list(obs), BAUs = BAUs, 
         nres = nres,
         response = "negative-binomial", link = "logit", 
         # manually set these arguments to reduce console output:
         K_type = "precision", method = "TMB", est_error = FALSE,
         normalise_wts = FALSE)

# ---- Predict over the BAUs ----

pred <- predict(S)

# ---- Plotting results ----


plot_list <- plot(S, pred$newdata, labels_from_coordnames = F)

plot_list <- lapply(plot_list, function(gg) {
    gg <- gg + xbreaks + ybreaks
    gg <- gg %>% change_font_size
    gg <- gg %>% change_legend_width
    gg <- gg %>% set_title_from_fill_legend
    return(gg)
  }
)

# Edit the legend breaks
if (simulation_method == "trig_field") {

  plot_list$p_prob  <- plot_list$p_prob %>% change_legend_breaks(breaks_prob)
  plot_list$interval90_prob <- plot_list$interval90_prob %>% change_legend_breaks(c(0.1, 0.15, 0.2))
  plot_list$p_mu <- plot_list$p_mu %>% change_legend_breaks(breaks_mu_BAU)
  plot_list$interval90_mu <- plot_list$interval90_mu %>% change_legend_breaks(c(50, 100, 150))
  
  
} else if (simulation_method == "model") {
  
  plot_list$interval90_prob <- plot_list$interval90_prob %>% change_legend_breaks(c(0.1, 0.175, 0.25))
  plot_list$p_prob <- change_legend_limits(plot_list$p_prob,  prob_scale_lims)
  plot_list$p_prob <- change_legend_breaks(plot_list$p_prob,  prob_scale_breaks)
  
}


figure <- create_figure_one_row_of_plots(list(plot_list$p_prob, 
                                              plot_list$interval90_prob, 
                                              plot_list$p_mu, 
                                              plot_list$interval90_mu))




ggsave(figure,
       filename = "3_2_Negbinom_sim_BAU_predictions.png", device = "png", 
       width = 14, height = 5.2,
       path = "results") 


# ---- Assess prediction and UQ at BAU level ----

## Compare BAUs$mu with the predictions (i.e., compute RMSPE), and the coverage
## via the percentiles.
pred$newdata$mu_true <- mu_true

# Only consider out-of-sample locations
unobsidx <- FRK:::unobserved_BAUs(S) 
pred$newdata <- pred$newdata[unobsidx, ] 

RMSPE <- function(z,pred) sqrt(mean((z - pred)^2))

## Coverage and MAPE
diagnostics <- pred$newdata@data %>% 
  mutate(
    true_in_pred_interval90 = mu_percentile_5 <= mu_true & mu_true <= mu_percentile_95,
    absolute_percentage_error = abs(p_mu - mu_true)/mu_true
  ) %>%
  dplyr::summarise(
    coverage_90 = mean(true_in_pred_interval90),
    MAPE = mean(absolute_percentage_error) * 100, 
    RMSPE = RMSPE(mu_true, p_mu)  
  ) 

## CRPS
diagnostics <- cbind(
  diagnostics,
  CRPS = mean(scoringRules::crps_sample(y = mu_true, dat = pred$MC$mu_samples))
)

write.csv(diagnostics, file = "Figures/Negbinom_sim.csv", row.names = FALSE)

save_html_table(
  diagnostics,
  file ="Figures/3_2_Negbinom_sim.html", 
  caption = "Negative-binomial simulation study"
)

# ---- Predict over polygons ----

pred_over_polygons <- predict(S, newdata = newdata)

plot_list <- plot(S, pred_over_polygons$newdata, colour = "black", labels_from_coordnames = F)

plot_list <- lapply(plot_list, function(gg) {
  gg <- gg + xbreaks + ybreaks
  gg <- gg %>% change_font_size
  gg <- gg %>% change_legend_width
  gg <- gg %>% set_title_from_fill_legend
  return(gg)
}
)

# Edit the legend breaks
if (simulation_method == "trig_field") {
  
  plot_list$p_prob  <- plot_list$p_prob %>% change_legend_breaks(c(0.5, 0.58, 0.66))
  plot_list$interval90_prob <- plot_list$interval90_prob %>% change_legend_breaks(c(0.016, 0.019, 0.022))
  plot_list$p_mu <- plot_list$p_mu %>% change_legend_breaks(c(8000, 11000, 14000))
  plot_list$interval90_mu <- plot_list$interval90_mu %>% change_legend_breaks(c(700, 1000, 1300))
  
  
} else if (simulation_method == "model") {
  
}


## Hack to make the correct x-axis and y-axis breaks appear
invisible_point <- geom_point(data = data.frame(x = c(0, 1), y = c(0, 1)), alpha = 0)
plot_list <- lapply(plot_list, function(gg) gg + invisible_point)

figure <- create_figure_one_row_of_plots(list(plot_list$p_prob, 
                                              plot_list$interval90_prob, 
                                              plot_list$p_mu, 
                                              plot_list$interval90_mu))

ggsave(figure,
       filename = "3_2_Negbinom_sim_arbitrary_polygon_predictions.png", device = "png", 
       width = 14, height = 5.2,
       path = "results") 

