## Define the functions we will use for comparing across methods
BrierScore <- function(Z, pred) mean((Z - pred)^2)
AUC <- function(Z, pred) as.numeric(pROC::auc(Z, pred))

## df should contain predicted probabilities in a field named pred, 
## and the validation data in a field called Z. 
compute_diagnostics_MODIS <- function(df) {
  summarise(df,
            Brier = BrierScore(Z, pred),
            AUC = AUC(Z, pred))
}


## Function to create long form dataframe, useful for diagnostics and plotting
long_prediction_data_frame <- function(df) {
  data.frame(
    Method = rep(PACKAGES, each = nrow(df)),
    x    = rep(df$x, times = length(PACKAGES)),
    y    = rep(df$y, times = length(PACKAGES)),
    pred   = c(as.matrix(df[, paste0("pred_", PACKAGES)]))
  )
}


## User defined functions for the MODIS comparison study

## Define a function needed for missing-block sampling.
## Samples a rectangular Block of width w and height h from a rectangular spatial 
## region defined by the coordinates x and y. 
## It randomly samples a point within D which acts as a 
## lower-left vertex for the Block. The sample region for this vertex depends on 
## the size of the Block; we cannot have it too close to the right- or top-edge of 
## D, otherwise the Block will extend beyond the considered boundaries. 
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


remove_geom <- function(ggplot2_object, geom_type) {
  # Delete layers that match the requested type.
  layers <- lapply(ggplot2_object$layers, function(x) {
    if (class(x$geom)[1] == geom_type) {
      NULL
    } else {
      x
    }
  })
  # Delete the unwanted layers.
  layers <- layers[!sapply(layers, is.null)]
  ggplot2_object$layers <- layers
  ggplot2_object
}


# ## For computing scoring rules, we need the MC samples of the test locations. 
# ## For producing plots, we need the MC samples at ALL locations.
# ## First, define a function to bind matrices with different numbers of columns.
# ## (not an ideal approach, but will suffice for now. I need to change this if 
# ## we use any plots that use the samples at all location).
# ## This function rbind() matrices together; if the numbre of columns is different, 
# ## it excludes the columns of the matrix with more columns.
# .concat_matrices <- function(X, Y) {
#   n_X <- ncol(X)
#   n_Y <- ncol(Y)
#   if (n_X != n_Y) {
#     print("Unequal number of columns; removing columns from matrix with excess columns.")
#     print(paste("Number of columns in X:", n_X))
#     print(paste("Number of columns in Y:", n_Y))
#     
#     
#     ## Compute discrepancy 
#     diff <- abs(n_Y - n_X)
#     ## Remove samples from which matrix had more samples (columns)
#     if (n_Y > n_X) {
#       Y <- Y[, -(n_Y:(n_Y - diff + 1))]
#     } else {
#       X <- X[, -(n_X:(n_X - diff + 1))]
#     }
#   } else {
#     print(paste("Number of columns equal. Number of columns:", n_X))
#   }
#   
#   ## Concatenate the matrices together
#   return(rbind(X, Y))
# }

