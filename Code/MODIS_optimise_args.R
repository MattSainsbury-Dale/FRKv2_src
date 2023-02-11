# devtools::install_github("MattSainsbury-Dale/testarguments")
# library("testarguments")
# setwd("~/Dropbox/FRKv2_src/")
source("Code/MODIS_diagnostic_and_misc_fns.R")

## Load the model fitting and prediction functions
invisible(
  mapply(source, paste0("Code/MODIS_modelling_fns/", PACKAGES, ".R"))
)


# ---- Construct validation set ----

## We do this by splitting the training set into two; a training set, 
## and a validation set: We do NOT use the test set to choose function 
## arguments. NB: df_train is defined in MODIS_analysis.R (using seed = 1 and MAR 
## sampling scheme).
idx_train <- 1:nrow(df_train)
RNGversion("3.6.0")
set.seed(1)
idx_fit <- sample(idx_train, round(nrow(df_train)/2), replace = FALSE)
idx_val <- idx_train[-idx_fit]

df_fit <- df_train[idx_fit, ]
df_val <- df_train[idx_val, ]

df_val$Z <- df_val$z

# ---- FRK ----

fun <- function(df_train, df_test, nres) {
  M <- MODIS_FRK_fit(df_train, nres = nres)
  pred <- MODIS_FRK_pred(df_test, M)
  
  ## test_arguments() requires fun to return a matrix-like object with named columns
  return(data.frame(pred = pred))
}

FRK_scores <- test_arguments(
  fun, df_fit, df_val, compute_diagnostics_MODIS, arguments = list(nres = 1:3)
)

plot_diagnostics(FRK_scores) 


# ---- spNNGP ----

fun <- function(df_train, df_test, n.neighbours) {
  M <- MODIS_spNNGP_fit(df_train, n.neighbours = n.neighbours)
  pred <- MODIS_spNNGP_pred(df_test, M)
  return(data.frame(pred = pred))
}

spNNGP_scores <- test_arguments(
  fun, df_fit, df_val, compute_diagnostics_MODIS, 
  arguments = list(n.neighbours = 1:20)
)

plot_diagnostics(spNNGP_scores) 



# ---- mgcv ----

fun <- function(df_train, df_test, k) {
  M <- MODIS_mgcv_fit(df_train, k = k)
  pred <- MODIS_mgcv_pred(df_test, M)
  return(data.frame(pred = pred))
}

mgcv_scores <- test_arguments(
  fun, df_fit, df_val, compute_diagnostics_MODIS, 
  arguments = list(k = c(
    seq(100, 1000, by = 100),
    seq(1200, 2000, by = 200), 
    seq(2250, 3000, by = 250) 
  ))
)

plot_diagnostics(mgcv_scores)



# ---- INLA ----

fun <- function(df_train, df_test, max.edge.interior) {
  M <- MODIS_mgcv_fit(df_train = df_train, pred_locs = df_test, max.edge.interior = max.edge.interior)
  return(data.frame(pred = pred))
}

INLA_scores <- test_arguments(
  fun, df_fit, df_val, compute_diagnostics_MODIS, 
  ## Run 1:
  # arguments = list(max.edge.interior = c(8, 10, 12),
  #                  sigma0 = c(3, 6, 9),
  #                  range0 = c(0.5, 1, 2))
  ## Run 2:
  arguments = list(max.edge.interior = seq(4, 8, by = 0.5))
)

plot_diagnostics(INLA_scores)



# ---- spBayes ----

fun <- function(df_train, df_test, knots) {
  M <- MODIS_spBayes_fit(df_train, knots = knots)
  pred <- MODIS_spBayes_pred(df_test, M)
  return(data.frame(pred = pred))
}

spBayes_scores <- test_arguments(
  fun, df_fit, df_val, compute_diagnostics_MODIS, 
  arguments = list(knots = (5:20)^2)
)

plot_diagnostics(spBayes_scores)


