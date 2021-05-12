library("FRK")
library("sp")
library("dplyr")
library("reshape2")

Am_data <- read.csv("./intermediates/Am_data.csv")
BAUs <- readRDS("./intermediates/Am_BAUs.rds")
blocks <- readRDS("./intermediates/Am_blocks.rds")


# ---- FRK ----

## In FRK, covariates are supplied at the BAU level: Remove covariates from the 
## data as the BAUs and the data cannot have common fields
Am_data[, c("x1", "x2", "x3", "x4")] <- NULL

coordinates(Am_data) = ~ Easting + Northing

BAUs$fs <- 1     # scalar matrix for fine scale variation
Am_data$std <- 1 # set measurement error to small value to replicate lognormal kriging

# M <- FRK(f = Am ~ x1 + x2 + x3, data = list(Am_data),
M <- FRK(f = Am ~ -1 + x1 + x2 + x3 + x4, data = list(Am_data),
         response = "gaussian", link = "log",
         BAUs = BAUs, est_error = FALSE)


# ---- Predict over the blocks ----

## Now we predict over arbitrary polygons.

## Predict using FRK
RNGversion("3.6.0")
set.seed(1)
pred <- predict(M, newdata = blocks, percentiles = NULL)

FRK_results <- pred$newdata@data
FRK_results$area_sqrt <- sapply(blocks@polygons, function(x) sqrt(x@area))
FRK_results$Scheme <- as.numeric(blocks@data$Scheme)
FRK_results$Framework <- "FRK" 

FRK_results <- FRK_results %>%
  dplyr::select(c("p_mu", "RMSPE_mu", "area_sqrt", "Scheme", "Framework")) %>% 
  melt(id.vars = c("area_sqrt", "Scheme", "Framework"))

write.csv(FRK_results, 
          file = "./intermediates/Am_FRK.csv", 
          row.names = FALSE)




