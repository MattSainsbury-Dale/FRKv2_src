
# ---- Create training and test sets ----

## First compute observed indices.
## The way we sample the observed indices depends on whether the sampling scheme
## is missing at random, or "missing at block".
if(!exists("seed")) {
  seed <- 1
}
RNGversion("3.6.0")
set.seed(seed)

if(!exists("missing_fwk")) {
  missing_fwk <- "MAR"
}
if(missing_fwk == "MAR") {
  n <- 6000
  obs.id <- sample(1:nrow(df), n, replace = FALSE) 
} else if (missing_fwk == "block") {
  block  <- sampleBlock(df, w = 30, h = 30)
  obs.id <- which(!with(df, (x > block["xmin"] & x <= block["xmax"]) & (y > block["ymin"] & y <= block["ymax"])))
}

## Based on the observed indices, construct the unobserved indices, 
## training set, and test set.
unobs.id <- (1:nrow(df))[-obs.id] 
df_train <- df[obs.id, ]   
df_test  <- df[unobs.id, ]


# ---- Optimise arguments ----

## See optimise_args.R and ../argument_selection/

# ---- Run the models ----

times         <- list()         ## store fitting times
fitted_model_objects <- list()  ## store model objects (which may be of use later)


## NB: MODIS_INLA() also assigns the fitted values (i.e., predictions at
## training locations) of the INLA to the parent environment (the environment
## from which the function was called) in an object named "INLA_fitted_values".
if("INLA" %in% PACKAGES)
  times$INLA <- system.time(
    df_test$pred_INLA <- MODIS_INLA(pred_locs = df_test, df_train = df_train, max.edge.interior = ARGS$max.edge.interior)
  )

if("FRK" %in% PACKAGES)
  times$FRK <- system.time({
    fitted_model_objects$FRK <- MODIS_FRK_fit(df_train, nres = ARGS$nres)
    df_test$pred_FRK <- MODIS_FRK_pred(df_test, fitted_model_objects$FRK)
  })

if("mgcv" %in% PACKAGES)
  times$mgcv <- system.time({
    fitted_model_objects$mgcv <- MODIS_mgcv_fit(df_train, k = ARGS$k)
    df_test$pred_mgcv <- MODIS_mgcv_pred(df_test, fitted_model_objects$mgcv)
  })

if("spNNGP" %in% PACKAGES)
  times$spNNGP <- system.time({
    fitted_model_objects$spNNGP <- MODIS_spNNGP_fit(df_train, n.neighbours = ARGS$n.neighbours)
    df_test$pred_spNNGP <- MODIS_spNNGP_pred(df_test, fitted_model_objects$spNNGP)
  })

if("spBayes" %in% PACKAGES)
  times$spBayes <- system.time({
    fitted_model_objects$spBayes <- MODIS_spBayes_fit(df_train, knots_squared = ARGS$knots_squared)
    df_test$pred_spBayes <- MODIS_spBayes_pred(df_test, fitted_model_objects$spBayes)
  })

times <- sapply(times, function(x) unname(x["elapsed"]))
times <- times[PACKAGES] 


# ---- Fitted values ----

## Compute fitted values so we can produce prediction maps over D

## Sanity check:
# all(row.names(fitted_model_objects$spNNGP$X) == obs.id)
# all(row.names(fitted_model_objects$spBayes$X) == obs.id)

if("INLA" %in% PACKAGES)
  df_train$pred_INLA <- INLA_fitted_values # NB: this is computed within MODIS_INLA()

## FRK doesn't have easily accessed fitted values.
## Just predict over the training locations using the fitted SRE object.
if("FRK" %in% PACKAGES)
  df_train$pred_FRK <- predict(
    fitted_model_objects$FRK,
    fitted_model_objects$FRK@data[[1]]
  )$newdata$p_prob

if("mgcv" %in% PACKAGES)
  df_train$pred_mgcv <- fitted_model_objects$mgcv$fitted.values

if("spNNGP" %in% PACKAGES)
  df_train$pred_spNNGP <- rowMeans(plogis(fitted_model_objects$spNNGP$y.hat.samples))

if("spBayes" %in% PACKAGES)
  df_train$pred_spBayes <- fitted_values_spBayes(fitted_model_objects$spBayes)


# ---- Output ----

## The controlling script now uses df_train and df_test, which now also contains
## fitted values and predictions, to compute diagnostics and ROC curves (using 
## df_test) and prediction maps over D (using rbind(df_train, df_test))

