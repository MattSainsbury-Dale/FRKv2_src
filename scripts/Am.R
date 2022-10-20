suppressMessages({
library("FRK")
library("sp")
library("ggplot2")
library("ggpubr")
library("dplyr")
library("reshape2")
library("DHARMa")

# georob has been archived since it depends on RandomFields; now, we load 
# the pre-saved estimates from georob.
# library("georob") 
# This package is required by georob, and it doesn't always install if
# we don't explicitly add it to dependencies.txt.
# library("maps")
})

# ---- Load and pre-process the Americium data ----

GZ_df   <- data.frame("Easting" = 219868.09, "Northing" = 285320.84)
Am_data <- read.csv("data/Am_data.csv")

## Convert Easting and Northing from feet to metres, and rename Americium.
Am_data$Easting   <- Am_data$Easting * 0.3048
Am_data$Northing  <- Am_data$Northing * 0.3048
Am_data$Am        <- Am_data$Americium
Am_data$Level     <- Am_data$Americium <- NULL

## Add covariates required for plotting figures from Paul and Cressie (2011), 
## and for prediction with georob.
d_cutoff <- 30.48  
d <- FRK::distR(Am_data[, c("Easting", "Northing")], GZ_df)
Am_data$d <- d # for use in Paul and Cressie
Am_data$x1 <- as.numeric(d < d_cutoff)
Am_data$x2 <- d * Am_data$x1
Am_data$x3 <- as.numeric(d >= d_cutoff)
Am_data$x4 <- d * Am_data$x3


# ---- Construct the blocking schemes  ----

# centre, width, and height
makeRectangle <- function(centre, w, h) {
  vertices <- rbind(c(centre[, 1] - w/2, centre[, 2] - h/2),
                    c(centre[, 1] - w/2, centre[, 2] + h/2),
                    c(centre[, 1] + w/2, centre[, 2] + h/2),
                    c(centre[, 1] + w/2, centre[, 2] - h/2),
                    c(centre[, 1] - w/2, centre[, 2] - h/2))
  Polygon(vertices)
}


## Following Paul and Cressie (2011),
## we predict over a series of concentric square blocks centred at Ground Zero 
## (GZ), as well as a series of concentric square blocks away from GZ.
construct_block_scheme <- function(n_schemes = 2, n_block = 5) {

  ratio <- 1.03  # width to height ratio of the blocks
  w     <- seq(43, 250, length.out = n_block)
  h     <- w / ratio
  
  ## Treat GZ as the centre, and expand relative to GZ.
  blocks <- list()
  for(i in 1:n_block) {
    blocks[[i]] <- makeRectangle(centre = GZ_df, w = w[i], h = h[i])
    blocks[[i]] <- Polygons(list(blocks[[i]]), paste0("block", i))
  }
  
  if (n_schemes == 2) {
    ## Now shift away from GZ
    centre <- GZ_df
    centre[, 1] <- GZ_df[, 1] - 153
    centre[, 2] <- GZ_df[, 2] + 125
    for(i in (n_block + 1):(2 * n_block)) {
      blocks[[i]] <- makeRectangle(centre = centre, w = w[i - n_block], h = h[i- n_block])
      blocks[[i]] <- Polygons(list(blocks[[i]]), paste0("block", i))
    }
  }
  
  
  ## (set the plotting order from largest to smallest)
  pred_polygons <- SpatialPolygons(blocks, (n_schemes * n_block):1)
  coordnames(pred_polygons) <- c("Easting", "Northing")
  
  pred_polygons$Scheme <- rep(as.character(n_schemes:1), each = length(pred_polygons)/n_schemes) 
  
  return(pred_polygons)
}

blocks <- construct_block_scheme()



# ---- Plot the data and the blocks ---- 

nasa_palette <- c(
  "#03006d","#02008f","#0000b6","#0001ef","#0000f6","#0428f6","#0b53f7",
  "#0f81f3","#18b1f5","#1ff0f7","#27fada","#3efaa3","#5dfc7b","#85fd4e",
  "#aefc2a","#e9fc0d","#f6da0c","#f5a009","#f6780a","#f34a09","#f2210a",
  "#f50008","#d90009","#a80109","#730005"
  )

lab1 <- xlab(as.expression(bquote("Easting /" ~ 10^5 ~ "m")))
lab2 <- ylab(as.expression(bquote("Northing /" ~ 10^5 ~ "m")))

## Basic plot to reduce code repetition
p_basic <- ggplot(data = Am_data, 
                  aes(x = Easting / 10^5, y = Northing / 10^5)) +
  lab1 + lab2 + scale_x_continuous(breaks = c(2.197, 2.199, 2.201)) + 
  scale_y_continuous(breaks = c(2.852, 2.854, 2.856)) + theme_bw() + coord_fixed()

## Data on the original scale
p_data <- p_basic +
  geom_point(aes(colour = Am), size = 1)  +
  geom_point(data = GZ_df, shape = 4, size = 5) +
  scale_colour_gradientn(colours = nasa_palette,
                         name = "Americium", 
                         labels = scales::scientific, 
                         breaks = c(250000, 750000))

## Data on the log scale
p_data_log_scale <- p_basic +
  geom_point(aes(colour = log(Am)), size = 1) +
  geom_point(data = GZ_df, shape = 4, size = 5) +
  scale_colour_gradientn(colours = nasa_palette,
                         name = "Log-Americium", 
                         breaks = c(9, 11, 13))

## Blocking scheme
p_Scheme_1_2 <- p_basic +
  # include the points with complete-transparancy so that the domain remains 
  # the same as the other two plots
  geom_point(size = 0.3, alpha = 0) + 
  geom_point(data = GZ_df, shape = 4, size = 5) +
  geom_polygon(data = FRK::SpatialPolygonsDataFrame_to_df(blocks), 
               aes(group = id, lty = Scheme), colour = "black", alpha = 0) +
  labs(colour = "Blocking Scheme")

ggsave( 
  ggpubr::ggarrange(p_data +
                      theme(legend.text=element_text(angle = 20)) + 
                      theme(text = element_text(size=17)), 
                    p_data_log_scale + theme(text = element_text(size=17)), 
                    p_Scheme_1_2 + theme(text = element_text(size=17)), 
                    nrow = 1, align = "hv", legend = "top"),
  filename = "4_2_Am_data_and_blocks.png", device = "png", width = 13.6, height = 4.5,
  path = "results/"
)


# ---- Construct the BAUs used for both FRK and georob ----

Am_df <- Am_data # keep a copy of Am data as a data frame
coordinates(Am_data) = ~ Easting + Northing
BAUs <- FRK::auto_BAUs(manifold = plane(), type = "grid", 
                       data = Am_data, nonconvex_hull = FALSE)

## Add covariates to the BAUs
d_cutoff <- 30.48  
d_BAU <- distR(coordinates(BAUs), GZ_df)
BAUs$x1 <- as.numeric(d_BAU < d_cutoff)
BAUs$x2 <- d_BAU * BAUs$x1
BAUs$x3 <- as.numeric(d_BAU >= d_cutoff)
BAUs$x4 <- d_BAU * (BAUs$x3)


# ---- georob analysis ---- 

# NB: This code is based on vignette("georob_vignette") 

# cat("Starting georob analysis...\n")
# 
# # Fit a spatial linear model by Gaussian (RE)ML 
# r.georob.m0.spher.reml <- georob(
#   log(Am) ~ -1 + x1 + x2 + x3 + x4, 
#   data = Am_df, 
#   locations = ~ Easting + Northing,
#   variogram.model = "RMexp",
#   param = c(variance = 0.1, nugget = 0.05, scale = 100), tuning.psi = 1000
#   )
# 
# # The diagnostics at the beginning of the summary output suggest that 
# # maximization of the restricted log-likelihood by nlminb() was successful. 
# 
# ## In the vignette, they refit the model with maximum likelihood in order to 
# ## perform step-wise covariate selection. Although we do not wish to perform
# ## step-wise covariate selection, we will still refit the model for consistency.
# r.georob.m0.spher.ml <- update(r.georob.m0.spher.reml,
#                                control=control.georob(ml.method="ML"))
# 
# 
# ## Lognormal block Kriging
# 
# ## If newdata is a SpatialPolygonsDataFrame then predict.georob() computes block
# ## Kriging predictions.
# ## However, first we need the covariates. georob requires one covariate value 
# ## for each polygon. To deal with this, we will average the BAU-level covariates
# ## within each polygon: 
# poly <- lapply(blocks@polygons, function(x) SpatialPolygons(list(x)))
# ind <- lapply(poly, function(x) over(as(BAUs, "SpatialPoints"), x))
# blocks$x1 <- sapply(ind, function(y) tapply(BAUs$x1, y, mean))
# blocks$x2 <- sapply(ind, function(y) tapply(BAUs$x2, y, mean))
# blocks$x3 <- sapply(ind, function(y) tapply(BAUs$x3, y, mean))
# blocks$x4 <- sapply(ind, function(y) tapply(BAUs$x4, y, mean))
# 
# ## Permanence of log-normality, that is, the assumption that both point values 
# ## and block means follow log-normal laws, strictly cannot hold. This does not 
# ## much impair the efficiency of the back-transformation as long as the blocks 
# ## are small (Cressie, 2006; Hofer et al., 2013). However, for larger blocks, 
# ## such as those used in this example, one should use the optimal predictor 
# ## obtained by averaging back-transformed point predictions. This is the approach
# ## that we take. 
# ## First, we need the full covariance matrix of the point prediction errors: To 
# ## Hence, we compute the point predictions of log(Am) with the additional 
# ## control argument full.covmat=TRUE:
# point_predictions <- predict(
#   r.georob.m0.spher.reml, newdata = as.data.frame(BAUs),
#   control = control.predict.georob(extended.output = TRUE, full.covmat = TRUE)
#   )
# 
# ## Now we back-transform the predictions and average them separately for each block:
# ## index defining to which block the point predictions belong
# poly <- lapply(blocks@polygons, function(x) SpatialPolygons(list(x)))
# ind <- lapply(poly, function(x) over(geometry(BAUs), x))
# ## select point predictions in block and predict block average
# block_predictions <- lapply(ind, function(i, x) {
#   idx <- which(!is.na(i))
#   x$pred <- x$pred[idx, ]
#   x$mse.pred <- x$mse.pred[idx, idx]
#   x$var.pred <- x$var.pred[idx, idx]
#   x$cov.pred.target <- x$cov.pred.target[idx, idx]
#   x$var.target <- x$var.target[idx, idx]
#   res <- lgnpp(x, is.block = TRUE)
#   return(res)
# }, x = point_predictions)
# block_predictions <- do.call(rbind, block_predictions)
# colnames(block_predictions) <- c("opt.pred", "opt.se")
# 
# ## Make a dataframe for use in the comparison plot
# georob_results <- block_predictions %>% 
#   as.data.frame() %>% 
#   rename(p_mu = opt.pred, RMSPE_mu = opt.se) %>% 
#   mutate(area_sqrt = sqrt(sapply(blocks@polygons, slot, "area")), 
#          Framework = "georob", 
#          Scheme = as.numeric(blocks@data$Scheme)) %>% 
#   melt(id.vars = c("area_sqrt", "Framework", "Scheme"))
# write.csv(as.data.frame(georob_results), file = "data/Am_georob_results.csv", row.names = FALSE)
# cat("georob analysis complete.\n")

# georob has been removed from CRAN since it depends on RandomFields
georob_results <- read.csv("data/Am_georob_results.csv") 


# ---- FRK analysis ---- 

cat("Starting FRK analysis...\n")

## In FRK, covariates are supplied at the BAU level: Remove covariates from the 
## data as the BAUs and the data cannot have common fields
Am_data@data[, c("x1", "x2", "x3", "x4")] <- NULL

# coordinates(Am_data) = ~ Easting + Northing

BAUs$fs <- 1     # scalar matrix for fine scale variation
Am_data$std <- 1 # set measurement error to small value to replicate lognormal kriging

suppressWarnings( # Suppress warning about log-link being inappropriate for Gaussian data 
  M <- FRK(f = Am ~ -1 + x1 + x2 + x3 + x4, data = list(Am_data),
           response = "gaussian", link = "log", BAUs = BAUs, est_error = FALSE, nres = 2,
           # manually set these arguments to reduce console output:
           K_type = "precision", method = "TMB")
)


## Predict over the blocks 
RNGversion("3.6.0"); set.seed(1)
pred <- predict(M, newdata = blocks, percentiles = NULL)

FRK_results <- pred$newdata@data %>% 
  mutate(area_sqrt = sapply(blocks@polygons, function(x) sqrt(x@area)), 
         Scheme = as.numeric(blocks@data$Scheme), 
         Framework = "FRK") %>%
  dplyr::select(c("p_mu", "RMSPE_mu", "area_sqrt", "Scheme", "Framework")) %>% 
  melt(id.vars = c("area_sqrt", "Scheme", "Framework"))

cat("FRK analysis complete.\n")


# ---- FRK diagnostic plots ----

DHARMA_object <- createDHARMa(
  simulatedResponse = simulate(M, conditional_fs = FALSE),
  observedResponse = Am_data$Am,
  integerResponse = FALSE
)

plot(DHARMA_object, quantreg = F)
plotResiduals(DHARMA_object, quantreg = F)
plotResiduals(DHARMA_object, form = sqrt(Am_data$d))


# ---- Only FRK ----

df <- FRK_results %>% 
  dplyr::select(c("variable", "value", "area_sqrt", "Scheme")) %>% 
  # alter the labels to change the facet_wrap titles:
  mutate(
    variable = factor(
    variable, 
    labels = c("'Block prediction'", "'RMSPE from block prediction'")
    )) %>% 
  mutate(
    Scheme = as.character(Scheme)
  )


figure <- ggplot(data = df, aes(x = area_sqrt, lty = Scheme)) +
  geom_line(aes(y = value), size = 1) +
  facet_wrap(~variable, scales = "free", labeller = label_parsed) + 
  labs(x = "Block size (m)", y = "", lty = "Blocking Scheme") +
  theme_bw() + 
  scale_y_continuous(labels = scales::scientific) + 
  theme(text = element_text(size = 20), 
        strip.text = element_text(size = 20))

figure <- figure + 
  theme(
    strip.background = element_blank()
  )

ggsave( 
  figure,
  filename = "4_2_Am_FRK.png", device = "png", width = 13.6, height = 4.5,
  path = "results/"
)



# ---- Comparison plot: FRK and georob ----


combined_df <- FRK_results %>% 
  dplyr::select(c("variable", "value", "area_sqrt", "Framework", "Scheme")) %>% 
  rbind(georob_results) %>% 
  mutate(fwk_sch = paste(Framework, Scheme, sep = ": ")) %>%
  # alter the labels to change the facet_wrap titles:
  mutate(variable = factor(
    variable, 
    labels = c("'Block prediction'", "'RMSPE from block prediction'")
  )) %>% 
  mutate(
    Framework = as.character(Framework), 
    Scheme = as.character(Scheme)
  )


figure <- ggplot(data = combined_df,
                 aes(x = area_sqrt, colour = Framework, 
                     lty = Scheme, group = fwk_sch)) +
  geom_line(aes(y = value), size = 1) +
  facet_wrap(~variable, scales = "free", labeller = label_parsed) + 
  labs(x = "Block size (m)", y = "", lty = "Blocking Scheme") +
  theme_bw() + 
  scale_y_continuous(labels = scales::scientific) + 
  theme(text = element_text(size = 20), 
        strip.text = element_text(size = 20))

ggsave( 
  figure,
  filename = "4_2_Am_comparison.png", device = "png", width = 13.6, height = 4.5,
  path = "results/"
)
