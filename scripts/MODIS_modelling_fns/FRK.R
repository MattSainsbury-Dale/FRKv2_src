MODIS_FRK_fit <- function(df_train, nres = 4) {
  
  ## Convert training data to SpatialPointsDataFrame
  zspdf <- df_train
  coordinates(zspdf) = ~ x + y     # Turn into sp object
  zspdf$k_Z <- 1 # constant size parameter of 1 (i.e., a Bernoulli distribution)  
  
  ## Construct Basis functions
  B <- auto_basis(data = zspdf, nres = nres)
  
  cat("FRK using", nres, "resolution(s) of basis functions,", nbasis(B), "basis functions in total.\n")
  
  ## BAUs (just use a grid over the spatial domain of interest)
  BAUs    <- SpatialPixelsDataFrame(points = expand.grid(x = 1:225, y = 1:150), 
                                    data = expand.grid(x = 1:225, y = 1:150))
  BAUs$fs <- 1   # Scalar variance matrix for fine-scale variation
  
  ## Construct SRE object
  M       <- SRE(f = z ~ 1, data = list(zspdf), basis = B, BAUs = BAUs, 
                 K_type = "precision", response = "binomial", link = "logit", 
                 # manually set these arguments to reduce console output:
                 est_error = FALSE)
  
  ## Fitting
  M <- SRE.fit(M, method = "TMB") 
  
  return(M)
}

MODIS_FRK_pred <- function(pred_locs, FRK_object) {
  
  ## Convert prediction locations to SpatialPointsDataFrame
  coordinates(pred_locs) = ~ x + y
  invisible(capture.output(
  pred <- predict(FRK_object, type = "mean", 
                  newdata = pred_locs,
                  percentiles = NULL)  
  ))
  
  return(pred$newdata$p_mu)
}
