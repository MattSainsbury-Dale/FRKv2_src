MODIS_mgcv_fit <- function(df_train, k = 2250) { 
  
  cat("mgcv using", k, "knots.\n")

  f <- z ~ te(x, y,            # inputs over which to smooth
              bs = "tp",       # type of bases
              k = k,           # knot count in each dimension
              d = 2)           # spatial basis dimension

  return(bam(f, family = binomial(link = "logit"), data = df_train))
}

MODIS_mgcv_pred <- function(pred_locs, mgcv_object) {
  
  return(predict(mgcv_object, newdata = pred_locs, type = "response"))
}



