# Here, we generate 1000 data points for each day of the year, so that we have
# 365,000 data points. We use poisson data. # TODO perhaps it'd be better to use inverse-gaussian or gamma data.
library("GpGp")
library("FRK")
library("spacetime")
library("sp")
library("ggplot2")
library("ggpubr")

RNGversion("3.6.0")
set.seed(1)
opts_FRK$set("verbose", TRUE)

# Spatio-temporal:
n <- 1000 * 365

lonlattime <-cbind(
  lon  = runif(n, min = -180, max = 180), 
  lat  = runif(n, min = -90, max = 90), 
  time = round(runif(n, min = 0, max = 365))
)

y  <- fast_Gp_sim(c(1, 20, 31, 0.5, 0), "matern_spheretime", lonlattime, m = 10)
y  <- 10 + y
mu <- y^2
Z <- rpois(length(mu), lambda = mu)

data <- STIDF(
  sp   = SpatialPoints(lonlattime[, c("lon", "lat")]),
  time = as.Date(lonlattime[, "time"], origin = "2020-01-01"),
  data = data.frame(Z = Z)
)

proj4string(data) <- CRS('+proj=longlat +ellps=sphere')


# ---- FRK analysis ----

# Automatically generate BAUs and basis functions
BAUs <- auto_BAUs(manifold = STsphere(), data = data, tunit = "days", isea3h_res = 5) 
BAUs$fs <- 1 # scalar fine-scale variance matrix
basis <- auto_basis(STplane(), data, tunit = "days", nres = 1) 


# Plot the spatial BAUs and the basis functions
SpatialBAUPlot <- plot_spatial_or_ST(
  BAUs@sp, column_names = "lon", alpha = 0, colour = "gray", 
  plot_over_world = TRUE)[[1]] + 
  theme(legend.position = "none")
SpatialBasisplot  <- show_basis(basis@Basis1)  
TemporalBasisplot <- show_basis(basis@Basis2) 

# Set up the model
S <- SRE(f = Z ~ 1, data = data, basis = basis, BAUs = BAUs,
         response = "poisson", link = "sqrt", K_type = "precision")

# Fit the model
S <- SRE.fit(S, method = "TMB")
saveRDS(S, file = "intermediates/Big_ST.rds")

# Predict over the BAUs
pred <- predict(S)

# Plotting
# l <- plot(S, pred, plot_over_world = TRUE, subset_time = 1:8) # TODO fix the bug that occurs when I run this; need to catch this error as it is so common! 
pred_time <- 367 # predict the next day
l <- plot(S, pred$newdata, plot_over_world = TRUE, subset_time = pred_time) 

# Apply a labeller so the facet shows day x rather than just x
facet_names <- paste0("Day ", pred_time)
names(facet_names) <- pred_time
l <- lapply(l, function(gg) gg + facet_wrap(~t, labeller = as_labeller(facet_names), nrow = 2))

gg <- ggarrange(l$p_mu, l$interval90_mu, align = "hv")
ggsave("Big_ST.png", plot = gg, path = "results", device = "png", width = 12, height = 3, units = "in") 

