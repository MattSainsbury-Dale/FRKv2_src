change_font_size <- function(gg) {
  gg + theme(axis.text = element_text(size = 16),
             axis.title = element_text(size = 19), 
             legend.text = element_text(size = 16), 
             plot.title = element_text(hjust = 0.5, size = 19))
}

change_legend_width <- function(gg, width = 1.1) {
  gg + theme(legend.key.width = unit(width, 'cm'))
}

## Remove y-axis labels/ticks for all but the left-most panel
interior_plot <- function(gg) {
  gg + rremove("ylab") + rremove("xlab") + rremove("y.text") + rremove("y.ticks")
}
exterior_plot <- function(gg) {
  gg + rremove("ylab") + rremove("xlab")
}

create_figure_one_row_of_plots <- function(plot_list) {
  plot_list[[1]] = plot_list[[1]] %>% exterior_plot
  for (i in 2:length(plot_list)) {
    plot_list[[i]] = plot_list[[i]] %>% interior_plot
  }
  
  figure <- ggarrange(plotlist = plot_list, 
                      nrow = 1, legend = "top", align = "hv") %>% 
    annotate_figure(left   = text_grob(bquote(s[2]), size = 20, vjust = 1,  hjust = 2, rot = 90),
                    bottom = text_grob(bquote(s[1]), size = 20, vjust = -1, hjust = -0.5))
} 


