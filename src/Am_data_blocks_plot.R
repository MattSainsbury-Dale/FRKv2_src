library("FRK")
library("ggplot2")
GZ_df <- read.csv("./intermediates/Am_GZ.csv")
Am_data <- read.csv("./intermediates/Am_data.csv")
blocks <- readRDS("./intermediates/Am_blocks.rds")



nasa_palette <- c("#03006d","#02008f","#0000b6","#0001ef","#0000f6","#0428f6","#0b53f7","#0f81f3",
                  "#18b1f5","#1ff0f7","#27fada","#3efaa3","#5dfc7b","#85fd4e","#aefc2a","#e9fc0d","#f6da0c","#f5a009",
                  "#f6780a","#f34a09","#f2210a","#f50008","#d90009","#a80109","#730005")

lab1 <- xlab(as.expression(bquote("Easting /" ~ 10^5 ~ "m")))
lab2 <- ylab(as.expression(bquote("Northing /" ~ 10^5 ~ "m")))

## Basic plot to reduce code repetition
p_basic <- ggplot(data = Am_data, 
                  aes(x = Easting / 10^5, y = Northing / 10^5)) +
  lab1 + lab2 + scale_x_continuous(breaks = c(2.197, 2.199, 2.201)) + 
  scale_y_continuous(breaks = c(2.852, 2.854, 2.856)) + theme_bw() + coord_fixed()

## Data on the original scale
p_data <- p_basic +
  geom_point(aes(colour = Am), size = 1)  +
  geom_point(data = GZ_df, shape = 4, size = 5) +
  scale_colour_gradientn(colours = nasa_palette,
                         name = "Americium", 
                         labels = scales::scientific, 
                         breaks = c(250000, 750000))

## Data on the log scale
p_data_log_scale <- p_basic +
  geom_point(aes(colour = log(Am)), size = 1) +
  geom_point(data = GZ_df, shape = 4, size = 5) +
  scale_colour_gradientn(colours = nasa_palette,
                         name = "Log-Americium", 
                         breaks = c(9, 11, 13))

## Blocking scheme
p_Scheme_1_2 <- p_basic +
  # include the points with complete-transparancy so that the domain remains 
  # the same as the other two plots
  geom_point(size = 0.3, alpha = 0) + 
  geom_point(data = GZ_df, shape = 4, size = 5) +
  geom_polygon(data = FRK::SpatialPolygonsDataFrame_to_df(blocks), 
               aes(group = id, colour = Scheme), alpha = 0) +
  labs(colour = "Blocking Scheme")

ggsave( 
  ggpubr::ggarrange(p_data + theme(legend.text=element_text(angle = 20)) + 
                      theme(text = element_text(size=17)), 
                    p_data_log_scale+ 
                      theme(text = element_text(size=17)), 
                    p_Scheme_1_2+ 
                      theme(text = element_text(size=17)), 
                    nrow = 1, align = "hv", legend = "top"),
  filename = "Am_data_and_blocks.png", device = "png", width = 13.6, height = 4.5,
  path = "./img/"
)
