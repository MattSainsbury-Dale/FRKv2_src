library("FRK")
library("sp")
library("ggplot2")
library("scoringRules") # crps_sample() 
library("dplyr")
library("ggpubr")

load("./data/Heaton_AllSatelliteTemps.RData")

## Create an identifier variable
df <- all.sat.temps %>%
  mutate(idx = 1:nrow(.)) 

## Find the id of observed and unobserved pixels
missing_idx <- filter(df, is.na(MaskTemp)) %>% pull(idx)
obs_idx     <- setdiff(1:nrow(df), missing_idx)
n_obs       <- length(obs_idx)

## Construct BAUs as SpatialPixels
BAUs <- SpatialPixelsDataFrame(points = df[, c("Lon", "Lat")], 
                               data = df[, c("Lon", "Lat")])
BAUs$fs <- 1   ## Fine-scale variation is iid

## Make training data as SpatialPoints
dat <- subset(df, !is.na(MaskTemp))   # No missing data in data frame
coordinates(dat)  <- ~Lon+Lat         # Convert to SpatialPointsDataFrame
dat$TrueTemp <- NULL                  # Remove TrueTemp

runtime <- system.time({
  ## Construct the basis functions
  basis <- auto_basis(plane(),            # we are on the plane
                      data = dat,         # data around which to make basis
                      regular = 1,        # regular basis
                      nres = 2,           # 3 resolutions
                      scale_aperture = 1) # aperture scaling of basis functions 
  
  ## Remove basis functions in problematic region
  if(nrow(df) == 150000) {
    basis_df <- data.frame(basis)
    rmidx <- which(basis_df$loc2 > 36.5 &
                     basis_df$loc1 > -94.5 &
                     basis_df$res == 3)
    basis <- remove_basis(basis, rmidx)
  }
  
  ## Construct SRE object
  M <- SRE(f = MaskTemp ~ 1 + Lon + Lat, data = list(dat), 
           basis = basis, BAUs = BAUs, K_type = "precision")
  
  ## Model fitting
  M <- SRE.fit(M, method = "TMB") 

  ## Prediction
  RNGversion("3.6.0")
  set.seed(1)
  pred <- predict(M, type = "response", percentiles = c(2.5, 97.5))
})

## Extract the dataframe from the Spatial object
pred_df <- pred$newdata@data

## Sanity Check: prediction and validation dataframes in same order
all(pred_df$Lon == df$Lon)
all(pred_df$Lat == df$Lat)

pred_df$TrueTemp <- df$TrueTemp
pred_df$id <- 1:nrow(pred_df)
validx <- which(!(pred_df$id %in% obs_idx) & !is.na(pred_df$TrueTemp))

## Sanity check:
nrow(na.omit(pred_df[-obs_idx, ])) == length(validx)

intervalScore <- function(Z, l, u, a = 0.05) {
  (u - l) + (2 / a) * (l - Z) * (Z < l) + (2 / a) * (Z - u) * (Z > u) 
}

## FIXME: Add run-time
write.csv(
  pred_df[validx, ] %>%
    summarise(RMSE = sqrt(mean((p_Z - TrueTemp)^2)),
              MAE = mean(abs(p_Z - TrueTemp)),
              CRPS = mean(crps_sample(y = pred_df[validx, "TrueTemp"], 
                                      dat = pred$MC$Z_samples[validx, ])),
              Cov95 = mean(Z_percentile_2.5 < TrueTemp & TrueTemp < Z_percentile_97.5), 
              intScore = mean(intervalScore(Z = TrueTemp, l = Z_percentile_2.5, u = Z_percentile_97.5))) %>% 
    as.data.frame(), 
          "./results/Heaton_FRKv2.csv", 
          row.names = FALSE)
