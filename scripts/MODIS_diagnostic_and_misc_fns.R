## Define the functions we will use for comparing across methods
BrierScore <- function(Z, pred) mean((Z - pred)^2)
AUC <- function(Z, pred) as.numeric(pROC::auc(Z, pred))

compute_diagnostics_MODIS <- function(df) {
  summarise(df, Brier = BrierScore(z, pred), AUC = AUC(z, pred))
}

## Creates long form dataframe, useful for diagnostics and plotting
long_prediction_data_frame <- function(df) {
  data.frame(
    Method = rep(PACKAGES, each = nrow(df)),
    x    = rep(df$x, times = length(PACKAGES)),
    y    = rep(df$y, times = length(PACKAGES)),
    pred   = c(as.matrix(df[, paste0("pred_", PACKAGES)]))
  )
}

## Define a function needed for the missing-block sampling scheme.
## Samples a rectangular Block of width w and height h from a rectangular spatial 
## region defined by the coordinates x and y in df. 
## It randomly samples a point within D which acts as a 
## lower-left vertex for the Block. The sample region for this vertex depends on 
## the size of the Block; we cannot have it too close to the right- or top-edge of 
## D, otherwise the Block will extend beyond the considered boundaries of D. 
sampleBlock <- function(df, w, h) {
  ## Retain only those coordinates which would result in the entire Block being inside D
  df <- df %>% subset(x < (max(x) - w) & y < (max(y) - h))
  
  ## Sample a point from the valid domain which will act as the lower-left vertex
  v1 <- df[sample(1:nrow(df), 1), c("x", "y")]
  
  ## Now create the upper-right vertex
  v2 <- v1 + c(w, h)
  
  return(c(xmin = v1[, 1], 
           ymin = v1[, 2], 
           xmax = v2[, 1], 
           ymax = v2[, 2]))
}

