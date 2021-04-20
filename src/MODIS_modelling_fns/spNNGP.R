## NB: If I can get n.omp.threads > 1 working, I can optimise the number of cores 
MODIS_spNNGP_fit <- function(df_train, n.neighbours = 15) {
  
  ## FIXME: 
  ## Strategy: the most important argument to optimise is n.neighbours. 
  ## Choose this first, then we can optimise n.samples and thin. 
  ## NB: thin is required at the prediction stage too... not sure how this will 
  ## work...
  ## Perhaps I can pass it through with the spNNGP object, if there is a slot 
  ## which is ignored in the current conditions. 
  ## can optimise this by optimising over n.samples
  # We want to use 400 samples in the end (for consistency with other packages)
  ## so solve: (n.samples - start)/thin = 400
  ## => start = n.samples - 400 * thin
  
  # ---- Fitting ----

  starting  <- list("phi" = 3/50, "sigma.sq" = 5)
  tuning    <- list("phi" = 0.15) ## Need to tweak this so that acceptance rate is between 25% and 50%
  priors    <- list("phi.Unif" = c(1/100, 10), "sigma.sq.IG" = c(2, 10)) # Changed this so prior for phi is more diffuse
  cov.model <- "exponential"
  
  ## samples to use
  n.samples <- 10000
  thin <- 10
  start <- n.samples - 400 * thin
  sub.sample <- list(start = start, thin = thin)
  
  M <- spNNGP::spNNGP(z ~ 1, coords = c("x", "y"), data = df_train, 
                          family = "binomial", method = "latent", 
                          n.neighbors = n.neighbours,
                          starting = starting, tuning = tuning, 
                          priors = priors, 
                          cov.model = cov.model, n.samples = n.samples, 
                          n.report = n.samples/2, 
                          sub.sample = sub.sample, 
                          fit.rep = TRUE, 
                          n.omp.threads = 1) # multi-core processing available through openMP

  return(M)
}
 
MODIS_spNNGP_pred <- function(pred_locs, spNNGP_object) {
  ## Predictions (out-of-sample testing locations only)
  ## could just pass sub-sample from the fitting function.
  n.samples <- nrow(spNNGP_object$p.beta.samples)
  thin <- 10
  start <- n.samples - 400 * thin
  sub.sample <- list(start = start, thin = thin)
  
  pred <- predict(spNNGP_object, 
                  X.0 = matrix(1, nrow = nrow(pred_locs)), 
                  coords.0 = cbind(x = pred_locs$x, y = pred_locs$y),
                  sub.sample = sub.sample,
                  n.report = n.samples/2)
  
  ## Predicted probabilities (testing locations)
  ## (note the logistic transformation)
  return(rowMeans(plogis(pred$p.y.0)))
}
