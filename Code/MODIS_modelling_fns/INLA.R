MODIS_INLA <- function(pred_locs, df_train, max.edge.interior = 5, 
                     sigma0 = 5, range0 = 1) {
  
  ## Establish a boundary for the domain, D
  boundary <- df_train[, c("x", "y")] %>% 
    slice_sample(n = min(2000, nrow(.))) %>%
    as.matrix() %>%
    inla.nonconvex.hull()
  
  ## Triangulation of the domain
  ## NB: Here I set the max edge length of the exterior of the domain to be 50%
  ## larger than the max interior length, and the cutoff to be 15% of the 
  ## maximum length. Reducing the number of arguments to check over in this way
  ## makes it easier to select a good configuration. 
  max.edge.exterior <- 1.5 * max.edge.interior
  cutoff <- 0.15 * max.edge.interior
  mesh <- inla.mesh.2d(boundary = boundary,
                       max.edge = c(max.edge.interior, max.edge.exterior), 
                       cutoff = cutoff)
  
  cat("INLA using a triangular mesh with", mesh$n, "vertices (and hence",  mesh$n,"basis functions)\n")
  
  ## Construct the SPDE on the mesh
  ## NB: I keep the probability associated with each prior (i.e., Psigma and 
  # Prange) fixed, and just optimise range0 and sigma0.
  spde <- inla.spde2.pcmatern(mesh = mesh, 
                              alpha = 2,
                              prior.range = c(range0, 0.01),
                              prior.sigma = c(sigma0, 0.01))
  
  ## Make an index which identifies spatial location 
  n_spatial <- mesh$n
  s_index   <- inla.spde.make.index(name = "spatial.field",
                                    n.spde = n_spatial)
  
  ## Fitting locations
  coords.fit <- df_train[ , c("x", "y")] %>% as.matrix()
  PHI <- inla.spde.make.A(mesh = mesh,
                          loc = coords.fit)
  
  ## First stack: Estimation
  n_data <- nrow(df_train)
  stack_est <- inla.stack(
    data    = list(z = df_train$z),
    A       = list(PHI, 1),
    effects = list(s_index,
                   list(Intercept = rep(1, n_data))),
    tag = "est")
  
  df_pred <- data.frame(x = mesh$loc[,1],
                        y = mesh$loc[,2])
  n_pred <- nrow(df_pred)
  PHI_pred <- Diagonal(n = n_pred)
  
  ## Second stack: Prediction
  stack_pred <- inla.stack(
    data = list(z = NA), # NA means we want to predict it
    A = list(PHI_pred, 1),
    effects = list(s_index,
                   list(Intercept = rep(1, n_pred))),
    tag = "pred")
  
  stack <- inla.stack(stack_est, stack_pred)
  
  ## Formula
  formula <- z ~ -1 + Intercept +
    f(spatial.field,
      model = spde,
      group = spatial.field.group)
  
  ## Fitting stage: almost all of the time is spent here
  inla.object <- inla(formula,
                     data = inla.stack.data(stack, spde = spde),
                     family = "binomial",
                     control.predictor = list(A = inla.stack.A(stack),
                                              compute = TRUE, 
                                              link = 1)) 
  
  ## Assign fitted values to the parent environment:
  index_est <- inla.stack.index(stack, "est")$data
  INLA_fitted_values    <- inla.object$summary.fitted.values$mean[index_est]
  assign("INLA_fitted_values", INLA_fitted_values, env = parent.frame()) 
  
  ## Prediction stage: this is comparatively fast
  ## Extract predictions and standard errors of basis function weights:
  index_pred <- inla.stack.index(stack, "pred")$data
  lp_mean    <- inla.object$summary.fitted.values$mean[index_pred]
  
  ## Locations we wish to predict over
  proj.grid <- inla.mesh.projector(
    mesh, loc = as.matrix(pred_locs[, c("x", "y")])
    )
  
  return(inla.mesh.project(proj.grid, lp_mean[1:mesh$n]))
}
