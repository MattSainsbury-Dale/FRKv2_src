library("georob")
library("sp")
library("dplyr")
library("reshape2")

Am_data <- read.csv("./intermediates/Am_data.csv")
BAUs <- readRDS("./intermediates/Am_BAUs.rds")
blocks <- readRDS("./intermediates/Am_blocks.rds")

# ---- 3.1 Exploratory analysis ----

## In the georob vignette, they do some exploratory analysis to find a 
## reasonable model. We will just use the model employed by Paul and Cressie 
## (2011): log(Americium) ~ -1 + x1 + x2 + x3 + x4

## Initial "drift" model
r.lm <- lm(log(Am) ~ -1 + x1 + x2 + x3 + x4, Am_data)

# ---- 3.2 Fitting a spatial linear model by Gaussian (RE)ML ----

# We fit the model that we developed in the exploratory analysis now by Gaussian REML:
r.georob.m0.spher.reml <- georob(log(Am) ~ -1 + x1 + x2 + x3 + x4, 
                                 Am_data, locations = ~Easting + Northing,
                                 variogram.model = "RMexp",
                                 param = c(variance = 0.1, nugget = 0.05, scale = 100),
                                 tuning.psi = 1000)

# The diagnostics at the begin of the summary output suggest that maximization of the
# restricted log-likelihood by nlminb() was successful. Nevertheless, before we interprete
# the output, we compute the profile log-likelihood for the range to see whether the maximization has found the global maximum:
r.prfl.m0.spher.reml.scale <- profilelogLik(r.georob.m0.spher.reml,
                                            values=data.frame(scale=seq(50, 200, by=10)))

## In the vignette, they refit the model with maximum likelihood in order to 
## perform step-wise covariate selection. Although we do not wish to perform
## step-wise covariate selection, we will still refit the model for consistency.
r.georob.m0.spher.ml <- update(r.georob.m0.spher.reml,
                               control=control.georob(ml.method="ML"))

## They do this; not sure it is necessary to update, but do it anyway. 
## Efficiency doesn't matter because we are not comparing run time.
r.georob.m1.spher.reml <- update(r.georob.m0.spher.reml)


# ---- 3.3 Computing Kriging predictions ----
# ---- 3.3.1 Lognormal point Kriging ----

## This section deals with prediction over SpatialPixelsDataFrame objects.
## That is, prediction over the BAUs.
coordinates(Am_data) <- ~ Easting + Northing

## computing Kriging predictions of log(zinc) and transforming them 
## back to the original scale of the measurements:
r.pk <- predict(r.georob.m0.spher.reml, newdata=BAUs,
                control=control.predict.georob(extended.output=TRUE))
r.pk <- lgnpp(r.pk)

# spplot(r.pk, zcol="lgn.pred", main="prediction")

pred <- spplot(r.pk, zcol="lgn.pred", main="prediction")
se <- spplot(r.pk, zcol="lgn.se", main="se")
# plot(pred, position=c(0, 0, 1/2, 1), more=TRUE)
# plot(se, position=c(1/2, 0, 1, 1), more=TRUE)


# ---- 3.3.2 Lognormal block Kriging ----

# If newdata is a SpatialPolygonsDataFrame then predict.georob() computes block
# Kriging predictions.

## Need the covariates. georob allows only one covariate value assigned to each polygon. 
## To deal with this, we will average over the polygon. 
## compute spatial mean of covariates for blocks:
poly <- lapply(blocks@polygons, function(x) SpatialPolygons(list(x)))
ind <- lapply(poly, function(x) over(as(BAUs, "SpatialPoints"), x))
blocks$x1 <- sapply(ind, function(y) tapply(BAUs$x1, y, mean))
blocks$x2 <- sapply(ind, function(y) tapply(BAUs$x2, y, mean))
blocks$x3 <- sapply(ind, function(y) tapply(BAUs$x3, y, mean))
blocks$x4 <- sapply(ind, function(y) tapply(BAUs$x4, y, mean))

## Block predictions with georob (assuming that both point values and block
## means follow log-normal laws — which strictly cannot hold):
r.bk <- predict(r.georob.m1.spher.reml, newdata = blocks,
                control=control.predict.georob(extended.output=TRUE, pwidth=25, pheight=25, mmax=25))
## transform the predictions back to the original scale by approximately 
## unbiased back-transformation proposed by Cressie (2006)
# newdata is used to pass point support covariates to lgnpp(), from which the spatial
# covariances of the covariates, needed for the back-transformation, are computed,
r.bk <- lgnpp(r.bk, newdata=BAUs)

pred <- spplot(r.bk, zcol="lgn.pred", main="prediction")
se <- spplot(r.bk, zcol="lgn.se", main="se")
# plot(pred, position=c(0, 0, 1/2, 1), more=TRUE)
# plot(se, position=c(1/2, 0, 1, 1), more=TRUE)

# The assumption that both point values and block means follow log-normal laws — which
# strictly cannot hold — does not much impair the efficiency of the back-transformation
# as long as the blocks are small (Cressie, 2006; Hofer et al., 2013). However, for larger
# blocks, one should use the optimal predictor obtained by averaging back-transformed point
# predictions. lgnpp() allows to compute this as well. 

## Optimal block predictions:
# We first compute block-Kriging predictions of log(zinc) and back-transform them as
# before under the permanence of log-normality assumption:
r.blks <- predict(r.georob.m1.spher.reml, newdata=blocks,
                  control=control.predict.georob(extended.output=TRUE, pwidth=800, pheight=800))
r.blks <- lgnpp(r.blks, newdata=BAUs)

# Note that we set pwidth and pheight equal to the dimension of the blocks. This is best
# for rectangular blocks because the blocks are then represented by a single pixel, which
# reduces computations.

# Next we use lgnpp() for computing the optimal log-normal block predictions. As we
# need the full covariance matrix of the point prediction errors for this we must re-
# compute the points predictions of log(zinc) with the additional control argument
# full.covmat=TRUE:
t.pk <- predict(r.georob.m0.spher.reml, newdata=as.data.frame(BAUs),
                control=control.predict.georob(extended.output=TRUE, full.covmat=TRUE))


# The predictions are now stored in the list component pred, and list components mse.pred,
# var.pred, var.target and cov.pred.target store the covariance matrices of prediction
# errors, predictions, prediction targets and the covariances between predictions and targets.
# Then we back-transform the predictions and average them separately for the blocks:


## index defining to which block the point predictions belong
poly <- lapply(blocks@polygons, function(x) SpatialPolygons(list(x)))
ind <- lapply(poly, function(x) over(geometry(BAUs), x))
## select point predictions in block and predict block average
tmp <- lapply(ind, function(i, x) {
  idx <- which(!is.na(i))
  x$pred <- x$pred[idx, ]
  x$mse.pred <- x$mse.pred[idx, idx]
  x$var.pred <- x$var.pred[idx, idx]
  x$cov.pred.target <- x$cov.pred.target[idx, idx]
  x$var.target <- x$var.target[idx, idx]
  res <- lgnpp(x, is.block = TRUE)
  return(res)
}, x = t.pk)
tmp <- do.call(rbind, tmp)
colnames(tmp) <- c("opt.pred", "opt.se")


# ---- Output ----

## All predictions
r.blks <- cbind(r.blks, tmp)

## Make a dataframe for use in the comparison plot
georob_results <- r.blks@data[, c("lgn.pred", "opt.pred", "lgn.se", "opt.se")]
georob_results <- data.frame(
  p_mu = c(georob_results$lgn.pred, georob_results$opt.pred), 
  RMSPE_mu = c(georob_results$lgn.se, georob_results$opt.se), 
  area_sqrt = sqrt(sapply(blocks@polygons, slot, "area")), 
  Framework = rep(c("georob: permanence of log-normality", 
                    "georob: optimal block prediction"), each = nrow(tmp)), 
  Scheme = rep(as.numeric(blocks@data$Scheme), times = 2)
)

georob_results <- georob_results %>% 
  melt(id.vars = c("area_sqrt", "Framework", "Scheme"))

write.csv(georob_results, 
          file = "./intermediates/Am_georob.csv", 
          row.names = FALSE)
