suppressMessages({
library("FRK")
library("gstat")
library("sp")
library("ggplot2")
library("scoringRules") # crps_sample() 
library("dplyr")
library("ggpubr")
source("Code/Utility_fns.R")
options(dplyr.summarise.inform = FALSE) # Suppress summarise info
})

## Use very-low-dimensional representations of the models to establish that the code works? 
quick <- check_quick()
nres <- if (quick) 2 else 4

cat(paste("Heaton study: Using", nres, "basis-function resolutions.\n"))

load("data/Heaton_AllSatelliteTemps.RData")

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
                      nres = nres,        # basis-function resolutions
                      scale_aperture = 1) # aperture scaling of basis functions 
  
  ## Remove basis functions in problematic region
  if(nrow(df) == 150000) {
    basis_df <- data.frame(basis)
    rmidx <- which(basis_df$loc2 > 36.5 &
                     basis_df$loc1 > -94.5 &
                     basis_df$res == 3)
    suppressMessages(basis <- remove_basis(basis, rmidx))
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

## Sanity Check (unit test): prediction and validation dataframes in same order
# all(pred_df$Lon == df$Lon)
# all(pred_df$Lat == df$Lat)

pred_df$TrueTemp <- df$TrueTemp
pred_df$id <- 1:nrow(pred_df)
validx <- which(!(pred_df$id %in% obs_idx) & !is.na(pred_df$TrueTemp))

## Sanity check (unit test):
# nrow(na.omit(pred_df[-obs_idx, ])) == length(validx)

intervalScore <- function(Z, l, u, a = 0.05) {
  (u - l) + (2 / a) * (l - Z) * (Z < l) + (2 / a) * (Z - u) * (Z > u) 
}

diagnostics <- pred_df[validx, ] %>%
  summarise(RMSE = sqrt(mean((p_Z - TrueTemp)^2)),
            MAE = mean(abs(p_Z - TrueTemp)),
            CRPS = mean(crps_sample(y = pred_df[validx, "TrueTemp"], 
                                    dat = pred$MC$Z_samples[validx, ])),
            Cov95 = mean(Z_percentile_2.5 < TrueTemp & TrueTemp < Z_percentile_97.5), 
            intScore = mean(intervalScore(Z = TrueTemp, l = Z_percentile_2.5, u = Z_percentile_97.5)), 
            runtime_minutes = runtime["elapsed"] / 60 # elapsed runtime in minutes
            ) %>% 
  as.data.frame()


write.csv(diagnostics, 
          "results/3_3_Heaton_FRKv2.csv", 
          row.names = FALSE)

diagnostics <- read.csv("results/3_3_Heaton_FRKv2.csv")

rownames(diagnostics) <- "FRK v2"

save_html_table(
  diagnostics,
  file = "results/3_3_Heaton_FRKv2.html", 
  caption = "Heaton comparison study"
)



# ---- Short (6-page) format material ----

short_format <- FALSE

if (short_format) {
  
  pred_FRKv2_df <- pred_df # rename it for clarity
  
  # ---- FRK v1 results ----
  
  # Same as above, but using only three resolutions of spatial basis functions; 
  # also use method = "EM".   
  
  ## Construct the basis functions
  basis <- auto_basis(plane(),            # we are on the plane
                      data = dat,         # data around which to make basis
                      regular = 1,        # regular basis
                      nres = 3,           # basis-function resolutions
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
  M_v1 <- SRE(f = MaskTemp ~ 1 + Lon + Lat, data = list(dat), 
              basis = basis, BAUs = BAUs, K_type = "block-exponential")
  
  ## Model fitting
  M_v1 <- SRE.fit(M_v1, method = "EM") 
  
  ## Prediction
  RNGversion("3.6.0")
  set.seed(1)
  pred_v1 <- predict(M_v1, type = "response")
  
  ## Extract the dataframe from the Spatial object
  pred_FRKv1_df <- pred_v1@data
  
  # ---- Plotting ----
  
  # Make a three-panel figure: Truth, FRK v1, and FRK v2. Zoom in on an 
  # interesting region that shows small-scale variation that can be captured with 
  # FRK v2, but not with FRK v1. 

  # Merge data frames:
  # NB: MaskTemp is the variable we used for training
  plot_df <- df %>% 
    select(Lon, Lat, TrueTemp) %>% 
    mutate(facet_var = "Truth") %>% 
    rename(variable = TrueTemp) %>% 
    rbind(
      pred_FRKv1_df %>% 
        select(Lon, Lat, mu) %>% 
        mutate(facet_var = "FRK v1") %>% 
        rename(variable = mu)
    ) %>% 
    rbind(
      pred_FRKv2_df %>% 
        select(Lon, Lat, p_Z) %>% 
        mutate(facet_var = "FRK v2") %>% 
        rename(variable = p_Z)
    )
  
  
  ## Save the data frame in case it is needed in the future:
  write.csv(plot_df, "results/Heaton_FRKv1_FRKv2.csv")
  plot_df <- read.csv("results/Heaton_FRKv1_FRKv2.csv")
  
  
  # Plot all:
  ggplot(plot_df, aes(x = Lon, y = Lat, fill = variable)) + 
    geom_raster() + 
    facet_wrap(. ~ facet_var) + 
    scale_fill_distiller(palette = "Spectral") + 
    theme_bw() + 
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0))
  
  
  # Focus on a window: 
  LonMin = -95.25; LonMax = -93
  LatMin = 35.9; LatMax = 36.5
  
  figure <- plot_df %>% 
    filter(LonMin <= Lon & Lon <= LonMax & LatMin <= Lat & Lat <= LatMax) %>% 
    ggplot(aes(x = Lon, y = Lat, fill = variable)) + 
    geom_raster() + 
    facet_wrap(. ~ facet_var, labeller = labeller(
      facet_var = labels
    )) + 
    scale_fill_distiller(palette = "Spectral", name = "") + 
    theme_bw() + 
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) + 
    theme(strip.background = element_blank(), strip.text.x = element_blank())
  
  nuke_axis <- function(gg) {
    gg + 
      rremove("ylab") + rremove("y.text") + rremove("y.ticks") +
      rremove("xlab") + rremove("x.text") + rremove("x.ticks")
  }
  
  figure <- figure %>% nuke_axis
  
  ggsave( 
    figure,
    filename = "Heaton_FRKv1_FRKv2.png", 
    device = "png", path = "results/", 
    width = 8, height = 1.7
  )
  
  # ---- Heaton table summary ----
  
  MAE <- c(1.33, 1.22, 1.65, 2.08, 1.33, 1.21, 1.24, 1.41, 2.05, 1.10, 1.87, 1.29)
  RMSPE <- c(1.86, 1.68, 2.08, 2.5, 1.85, 1.64, 1.68, 1.8, 2.52, 1.53, 2.45, 1.79)
  CRPS <- c(1.17, 0.87, 1.17, 1.44, 0.94, 0.85, 0.87, 1.02, 1.85, 0.83, 1.32, 0.91)
  Cvg95 <- c(0.36, 0.96, 0.83, 0.89, 0.92, 0.95, 0.94, 0.86, 0.75, 0.97, 0.93, 0.93)
  diagnostics <- data.frame(MAE = MAE, RMSPE = RMSPE, CRPS = CRPS, Cvg95 = Cvg95)
  apply(diagnostics, 2, mean)
  
}