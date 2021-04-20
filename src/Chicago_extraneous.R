# ---- User defined functions ----

## This script contains code which were previously used in the Chicago analysis, 
## and I keep in case they are useful in the future.


## Add zero counts to data (this was used when I had fine BAUs rather than the 
## large community areas)
add_zero_counts <- function (BAUs, data, type) {
  
  if (is(BAUs, "SpatialPixelsDataFrame")) {
    zero_counts       <- as.data.frame(coordinates(BAUs))
    zero_counts$count <- 0
    coordinates(zero_counts) <- coordnames(BAUs)
    # coordinates(zero_counts) = ~ longitude + latitude
    data   <- rbind(data, zero_counts)
    
  } else if (is(BAUs, "STFDF")) {
    zero_counts       <- as.data.frame(coordinates(BAUs))
    zero_counts$year  <- BAUs@data$t + 2000   
    zero_counts$count <- 0
    
    ## Convert the zero counts to STIDF, rbind() to the original data, then create the SRE object.
    zero_counts <- within(zero_counts, {time = as.Date(paste(year, 06,01,sep="-"))}) 
    zero_counts <- stConstruct(x = zero_counts, 
                               space = c("longitude", "latitude"), 
                               time = "time", 
                               interval = TRUE) 
    data   <- rbind(data, zero_counts)
    
  } 
  
  return(data)
}


## Determine the approximate area of a BAU cell
BAU_area <- function (BAUs) {
  
  if(!(is(BAUs, "SpatialPixelsDataFrame"))) 
    stop("Requires SpatialPixelsDataFrame")
  
  
  BAU_width <- as.data.frame(BAUs@grid)
  point1 <- BAU_width$cellcentre.offset
  point2 <- point1 + BAU_width$cellsize
  
  BAU_diag_length <- geosphere::distHaversine(point1, point2)
  ## Convert from m to km
  BAU_diag_length <- BAU_diag_length / 1000
  BAU_area <- (BAU_diag_length)^2 / 2
  
  return(BAU_area)
}



add_spatial_covariate_to_ST_BAUs <- function(ST_BAUs, binned_data, cov_name) {
  
  
  if (length(binned_data) != nrow(ST_BAUs@sp))
    warning("Some BAUs did not have any covariate data. Setting the data for these BAUs to zero.")
  
  ## Some BAUs may not have any population data. 
  ## Removing these BAUs makes the prediction maps look patchy and weird, so instead we will set 
  ## the population of these BAUs to 0. 
  ST_BAUs_spatial_frame <- ST_BAUs@sp
  ST_BAUs_spatial_frame@data[, cov_name] <- 0
  x <- rownames(ST_BAUs_spatial_frame@data)
  idx <- match(binned_data@data$BAU_name, x)
  ST_BAUs_spatial_frame@data[idx, cov_name] <- binned_data@data[, cov_name]
  
  ## Now replicate over each time
  tmp <- as.data.frame(rep(ST_BAUs_spatial_frame@data[, cov_name], length(ST_BAUs@time)))
  colnames(tmp) <- cov_name
  
  ST_BAUs@data <- cbind(
    ST_BAUs@data,
    tmp
  )
  
  return(ST_BAUs)
}