## Ground zero
GZ_df <- data.frame("Easting" = 219868.09, "Northing" = 285320.84)

write.csv(GZ_df, 
          file = "./intermediates/Am_GZ.csv", 
          row.names = FALSE)
