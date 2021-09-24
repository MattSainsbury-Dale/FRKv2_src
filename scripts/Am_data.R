suppressMessages(
library("FRK")
)
Am_data <- read.csv("./data/Am_data.csv")
GZ_df <- read.csv("./intermediates/Am_GZ.csv")

## Convert Easting and Northing from feet to metres, and rename Americium.
Am_data$Easting   <- Am_data$Easting * 0.3048
Am_data$Northing  <- Am_data$Northing * 0.3048
Am_data$Am        <- Am_data$Americium
Am_data$Level <- Am_data$Americium <- NULL

# save(Am_data, file="~/Dropbox/FRK/data/Am_data.rda")

## Add covariates required for plotting figures from Paul and Cressie (2011), 
## and for prediction with georob.
d_cutoff <- 30.48  
d <- FRK::distR(Am_data[, c("Easting", "Northing")], GZ_df)
Am_data$d <- d # for use in Paul and Cressie
Am_data$x1 <- as.numeric(d < d_cutoff)
Am_data$x2 <- d * Am_data$x1
Am_data$x3 <- as.numeric(d >= d_cutoff)
Am_data$x4 <- d * Am_data$x3

write.csv(Am_data, 
          file = "./intermediates/Am_data.csv", 
          row.names = FALSE)


