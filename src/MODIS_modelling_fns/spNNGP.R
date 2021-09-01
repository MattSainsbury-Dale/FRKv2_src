## NB: If I can get n.omp.threads > 1 working, I can optimise the number of cores 
MODIS_spNNGP_fit <- function(df_train, n.neighbours = 15) {
  
  # ---- Fitting ----

  starting  <- list("phi" = 3/50, "sigma.sq" = 5)
  tuning    <- list("phi" = 0.15) 
  priors    <- list("phi.Unif" = c(1/100, 10), "sigma.sq.IG" = c(2, 10))
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
                          n.omp.threads = 1)

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
