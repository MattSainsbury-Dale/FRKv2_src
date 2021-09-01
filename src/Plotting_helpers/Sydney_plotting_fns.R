## These plotting functions are saved here to simplify the main analysis script.

## Simplify the legend of the predictions if the fitting method was SA2s only 
## (I do this for very short 15 minute presentations where I do not have time 
## to explain everything properly and I don't want to overload the viewer) 
simplify_legend_label <- function(plots) {
    plots$p_prob <- plots$p_prob + 
      labs(fill = "Predicted proportion\nof families in poverty")
    plots$interval90_prob <- plots$interval90_prob + 
      labs(fill = "90% prediction-interval \nwidth for proportion of\nfamilies in poverty")
  
  return(plots)
}


change_font_size_and_axis <- function(gg) {
  gg + theme(axis.text = element_text(size = 11),
             axis.title = element_text(size = 14), 
             legend.text = element_text(size = 11),
             legend.title = element_text(size = 14)) + 
    scale_x_continuous(breaks = c(150.8, 151.0, 151.2), expand = c(0, 0)) + 
    scale_y_continuous(breaks = c(-33.7, -33.9, -34.1), expand = c(0, 0))
}