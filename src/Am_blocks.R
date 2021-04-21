library("sp")
GZ_df <- read.csv("./intermediates/Am_GZ.csv")

# centre, width, and height
makeRectangle <- function(centre, w, h) {
  vertices <- rbind(c(centre[, 1] - w/2, centre[, 2] - h/2),
                    c(centre[, 1] - w/2, centre[, 2] + h/2),
                    c(centre[, 1] + w/2, centre[, 2] + h/2),
                    c(centre[, 1] + w/2, centre[, 2] - h/2),
                    c(centre[, 1] - w/2, centre[, 2] - h/2))
  Polygon(vertices)
}


## Following Paul and Cressie (2011),
## we predict over a series of concentric square blocks centred at Ground Zero 
## (GZ), as well as a series of concentric square blocks away from GZ.
construct_block_scheme <- function() {
  n_schemes <- 2
  n_block <- 5
  
  ratio <- 1.03  # width to height ratio of the blocks
  w     <- seq(43, 250, length.out = n_block)
  h     <- w / ratio
  
  ## Treat GZ as the centre, and expand relative to GZ.
  blocks <- list()
  for(i in 1:n_block) {
    blocks[[i]] <- makeRectangle(centre = GZ_df, w = w[i], h = h[i])
    blocks[[i]] <- Polygons(list(blocks[[i]]), paste0("block", i))
  }
  
  ## Now shift away from GZ
  centre <- GZ_df
  centre[, 1] <- GZ_df[, 1] - 153
  centre[, 2] <- GZ_df[, 2] + 125
  for(i in (n_block + 1):(2 * n_block)) {
    blocks[[i]] <- makeRectangle(centre = centre, w = w[i - n_block], h = h[i- n_block])
    blocks[[i]] <- Polygons(list(blocks[[i]]), paste0("block", i))
  }
  
  ## (set the plotting order from largest to smallest)
  pred_polygons <- SpatialPolygons(blocks, (2 * n_block):1)
  coordnames(pred_polygons) <- c("Easting", "Northing")
  
  pred_polygons$Scheme <- rep(c("2", "1"), each = length(pred_polygons)/2) 
  
  return(pred_polygons)
}

blocks <- construct_block_scheme()

saveRDS(blocks, "./intermediates/Am_blocks.rds")
