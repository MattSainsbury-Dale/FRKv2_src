change_legend_breaks <- function(gg, breaks, aesthetic = "fill") {
  
  ## Find the scales associated with the specifed aesthetic
  sc <- as.list(gg$scales)$scales
  all_aesthetics <- sapply(sc, function(x) x[["aesthetics"]][1]) 
  idx <- which(aesthetic == all_aesthetics) 
  
  ## Overwrite the breaks of the specifed aesthetic
  gg$scales$scales[[idx]][["breaks"]] <- breaks
  
  return(gg)
} 

# same as above but change the limits
change_legend_limits <- function(gg, limits, aesthetic = "fill") {
  
  ## Find the scales associated with the specifed aesthetic
  sc <- as.list(gg$scales)$scales
  all_aesthetics <- sapply(sc, function(x) x[["aesthetics"]][1]) 
  idx <- which(aesthetic == all_aesthetics) 
  
  ## Overwrite the breaks of the specifed aesthetic
  gg$scales$scales[[idx]][["limits"]] <- limits
  
  return(gg)
} 



set_title_from_fill_legend <- function(gg) {
  gg + labs(title = gg$labels$fill) + labs(fill = "")
}