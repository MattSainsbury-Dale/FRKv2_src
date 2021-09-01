cloud_colour    = "orange"  # "white"
no_cloud_colour = "blue"    # "black"
missing_colour  = "white"   # "#BFD5E3"
midpoint_colour = "gray" 

change_font_size <- function(gg) {
  gg + theme(axis.text = element_text(size = 10),
             axis.title = element_text(size = 13), 
             legend.text = element_text(size = 10), 
             strip.text = element_text(size = 12))
}

training_data_background <- theme(
  panel.background = element_rect(fill = "white", colour = "white"), # "#6D9EC1"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

discrete_cloud_theme <- theme(legend.key=element_rect(
  # fill = "gray", colour = "gray", # use this if you colour the clouds white
  fill = "white", colour = "white", 
  size = 2))


discrete_cloud_scale <- scale_fill_gradient(low = no_cloud_colour, 
                                            high = cloud_colour,
                                            breaks = c(0, 1),
                                            guide = "legend",
                                            labels = c("No Cloud", "Cloud"),
                                            name = "")


common_layers <- ggplot() + theme_bw() + coord_fixed() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + 
  labs(x = bquote(s[1]), y = bquote(s[2]))




