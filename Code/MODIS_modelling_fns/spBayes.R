## NB: unlike spNNGP, the returned fitted values for y are on the probability
## scale, and so we do not need to apply the logistic transformation.
MODIS_spBayes_fit <- function(df_train, knots = 20^2) {
  
  if(sqrt(knots) != round(sqrt(knots)))
    stop("Number of knots must be a square number")
  
  cat("spBayes using", knots, "knots.\n")
  
  
  # ---- Fitting ----
  
  ## Get starting parameter values and tuning parameters
  fit <- glm(z ~ 1, data = df_train, family = "binomial")
  beta.starting <- coefficients(fit)
  beta.tuning   <- t(chol(vcov(fit)))
  
  ## Parameters we can change to make the model faster
  n.batch      <- 200 
  batch.length <- 50  

  ## knots greatly impacts fitting time - 
  ## prod(knots) is the total number of knots
  ## (note that the spNNGP paper compared the performance of spNNGP against 
  ## spGLM with only 25 knots! I think that may be too low and not fair to spGLM).
  knots   <- sqrt(c(knots, knots))
  
  ## Fit the model
  M <- spBayes::spGLM(z ~ 1, family = "binomial", data = df_train, 
                        coords = as.matrix(df_train[, c("x", "y")]), 
                        starting = list("beta" = beta.starting, "phi" = 0.06, "sigma.sq" = 1, "w" = 0),
                        tuning = list("beta" = beta.tuning, "phi" = 0.5, "sigma.sq" = 0.5, "w" = 0.5),
                        priors = list("beta.Normal" = list(0,10), "phi.Unif" = c(0.03, 0.3), "sigma.sq.IG" = c(2, 1)),
                        amcmc = list("n.batch" = n.batch, "batch.length" = batch.length, "accept.rate" = 0.43),
                        cov.model = "exponential", verbose = TRUE, n.report = n.batch/4, 
                        knots = knots)  
  return(M)
}

MODIS_spBayes_pred <- function(pred_locs, spBayes_object) {
  
  n.samples <- nrow(spBayes_object$p.beta.theta.samples)
  burn.in   <- 0.6001 * n.samples 
  
  ## Predictions (out-of-sample testing locations only)
  pred <- spBayes::spPredict(spBayes_object, 
                             pred.coords = as.matrix(pred_locs[, c("x", "y")]), 
                             pred.covars = matrix(1, nrow = nrow(pred_locs)),
                             start = burn.in,
                             thin = 10,
                             n.report = 500)
  
  ## Predicted probabilities (testing locations)
  return(rowMeans(pred$p.y.predictive.samples))
}

fitted_values_spBayes <- function(M) {
  n.samples <- nrow(M$p.beta.theta.samples)
  burn.in   <- 0.6001 * n.samples
  sub.samps <- burn.in:n.samples
  thin <- 10
  sub.samps <- sub.samps[seq(1, length(sub.samps), thin)]
  
  beta.hat <- M$p.beta.theta.samples[sub.samps,"(Intercept)"]
  w.hat    <- M$p.w.samples[, sub.samps]
  
  ## fitted probabilities (training locations)
  p.hat <- 1 / (1 + exp(-(M$X %*% beta.hat + w.hat)))
  p.hat.mu <- rowMeans(p.hat)
  
  return(p.hat.mu)
}