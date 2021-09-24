suppressMessages({
library("FRK")
library("sp")
})
  
Am_data <- read.csv("./intermediates/Am_data.csv")
GZ_df <- read.csv("./intermediates/Am_GZ.csv")

coordinates(Am_data) = ~ Easting + Northing

BAUs <- FRK::auto_BAUs(manifold = plane(), type = "grid", 
                       data = Am_data, nonconvex_hull = FALSE)

## Add covariates to the BAUs
d_cutoff <- 30.48  
d_BAU <- distR(coordinates(BAUs), GZ_df)

## Using only three covariates (as per the manuscript; it leads to the same results, 
## because they are equivalent models, but I wrote the georob and Paul and Cressie scripts
## using x1 to x4, and I dont want to change it)
# BAUs$x1 <- d_BAU * (d_BAU < d_cutoff)
# BAUs$x2 <- d_BAU >= d_cutoff
# BAUs$x3 <- d_BAU * (d_BAU >= d_cutoff)

BAUs$x1 <- as.numeric(d_BAU < d_cutoff)
BAUs$x2 <- d_BAU * BAUs$x1
BAUs$x3 <- as.numeric(d_BAU >= d_cutoff)
BAUs$x4 <- d_BAU * (BAUs$x3)

saveRDS(BAUs, "./intermediates/Am_BAUs.rds")
