## This script creates df, the MODIS cloud dataframe used in the FRK v2 paper
## It first lowers the resolution of the image (from over 10 million pixels to 
## 33750), and then applies a threshold to convert from continuous to binary data.

## Load the full, original dataset.
## This dataset is called df, and contains over 10 million pixels.
## We do not need this many, so we will coarsen the dataset in this script.
load("~/Dropbox/FRK/Data/MODISFullData.rda")
nrow(df)  # original number of pixels

glimpse(df)

apply(df, 2, range)
unique(df$x)


## Construct coarse boxes (150 x 225 grid)
## Choose x to contain 1.5 times as many cutpoints because
## there are approx 1.5 times as many x values in original dataset.
df$xcut <- cut(df$x, breaks = 225, labels = FALSE)
df$ycut <- cut(df$y, breaks = 150, labels = FALSE)

## Compute the mean of 'z' for each coarse box and aggregate original data 
## into lower resolution pixels
df <- aggregate(z ~ xcut + ycut, data = df, FUN = "mean" )
df <- dplyr::rename(df, x = xcut)
df <- dplyr::rename(df, y = ycut)
## We now have 33750 pixels instead of over 10 million.

## Threshold the response variable (radiance -> cloud/no cloud)
## Apply a threshold on 'z' for cloud/no cloud.
## cloud = 1 if z > thresh, cloud = 0 if z <= thresh.
df$z_unthresholded <- df$z 
thresh <- 7000
df$z <- 1 * (df$z > thresh)
rm(thresh)

## IMPORTANT:
## The name of the data frame object must also be MODIS_cloud_df when we save it.
MODIS_cloud_df <- df
rm(df)
rm(.Random.seed)

## Use version argument to specify a workspace format that is supported for 
## the R version which FRK currently depends (3.2.0). The default version used
## by save() is version 3, which is only supported by R 3.5.0 and
## above. Using this version adds a dependency to R version >= 3.5.0.
## If FRK ever explicitly depends on R >= 3.5.0, we are free to use the defult 
## version 3 in these save functions. 
# save(MODIS_cloud_df, file="MODIS_cloud_df.rda", version = 2)
save(MODIS_cloud_df, file="MODIS_cloud_df.rda")

