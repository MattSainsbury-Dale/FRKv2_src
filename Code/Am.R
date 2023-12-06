suppressMessages({
library("FRK")
library("sp")
library("ggplot2")
library("ggpubr")
library("dplyr")
library("reshape2")
})

# ---- Load and pre-process the Americium data ----

GZ_df   <- data.frame("Easting" = 219868.09, "Northing" = 285320.84)
Am_data <- read.csv("data/Am_data.csv")
write.csv(sd(Am_data$Americium), file = "Figures/Am_total_data_standard-deviation.csv", row.names = FALSE)


## Convert Easting and Northing from feet to metres, and rename Americium.
Am_data$Easting   <- Am_data$Easting * 0.3048
Am_data$Northing  <- Am_data$Northing * 0.3048
Am_data$Am        <- Am_data$Americium
Am_data$Level     <- Am_data$Americium <- NULL

## Add covariates required for plotting figures from Paul and Cressie (2011), 
## and for prediction with georob.
d_cutoff <- 30.48  
d <- FRK::distR(Am_data[, c("Easting", "Northing")], GZ_df)
Am_data$d <- d # for use in Paul and Cressie
Am_data$x1 <- as.numeric(d < d_cutoff)
Am_data$x2 <- d * Am_data$x1
Am_data$x3 <- as.numeric(d >= d_cutoff)
Am_data$x4 <- d * Am_data$x3


# ---- Construct the blocking schemes  ----

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
construct_block_scheme <- function(n_schemes = 2, n_block = 5) {

  ratio <- 1.03  # width to height ratio of the blocks
  w     <- seq(43, 250, length.out = n_block)
  h     <- w / ratio
  
  ## Treat GZ as the centre, and expand relative to GZ.
  blocks <- list()
  for(i in 1:n_block) {
    blocks[[i]] <- makeRectangle(centre = GZ_df, w = w[i], h = h[i])
    blocks[[i]] <- Polygons(list(blocks[[i]]), paste0("block", i))
  }
  
  if (n_schemes == 2) {
    ## Now shift away from GZ
    centre <- GZ_df
    centre[, 1] <- GZ_df[, 1] - 153
    centre[, 2] <- GZ_df[, 2] + 125
    for(i in (n_block + 1):(2 * n_block)) {
      blocks[[i]] <- makeRectangle(centre = centre, w = w[i - n_block], h = h[i- n_block])
      blocks[[i]] <- Polygons(list(blocks[[i]]), paste0("block", i))
    }
  }
  
  
  ## (set the plotting order from largest to smallest)
  pred_polygons <- SpatialPolygons(blocks, (n_schemes * n_block):1)
  coordnames(pred_polygons) <- c("Easting", "Northing")
  
  pred_polygons$Scheme <- rep(as.character(n_schemes:1), each = length(pred_polygons)/n_schemes) 
  
  return(pred_polygons)
}

blocks <- construct_block_scheme(n_schemes = 1)



# ---- Plot the data and the blocks ---- 

nasa_palette <- c(
  "#03006d","#02008f","#0000b6","#0001ef","#0000f6","#0428f6","#0b53f7",
  "#0f81f3","#18b1f5","#1ff0f7","#27fada","#3efaa3","#5dfc7b","#85fd4e",
  "#aefc2a","#e9fc0d","#f6da0c","#f5a009","#f6780a","#f34a09","#f2210a",
  "#f50008","#d90009","#a80109","#730005"
  )

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
  geom_point(size = 0.3, alpha = 1) + 
  geom_point(data = GZ_df, shape = 4, size = 5) +
  geom_polygon(data = FRK::SpatialPolygonsDataFrame_to_df(blocks), 
               aes(group = id, lty = Scheme), colour = "black", alpha = 0) +
  labs(lty = "Blocks")

ggsave( 
  ggpubr::ggarrange(p_data +
                      theme(legend.text=element_text(angle = 20)) + 
                      theme(text = element_text(size=17)), 
                    p_data_log_scale + theme(text = element_text(size=17)), 
                    p_Scheme_1_2 + theme(text = element_text(size=17)), 
                    nrow = 1, align = "hv", legend = "top"),
  filename = "4_2_Am_data_and_blocks.png", device = "png", width = 13.6, height = 4.5,
  path = "Figures/"
)


# ---- Construct the BAUs used for both FRK and georob ----

Am_df <- Am_data # keep a copy of Am data as a data frame
coordinates(Am_data) = ~ Easting + Northing
BAUs <- FRK::auto_BAUs(manifold = plane(), type = "grid", 
                       data = Am_data, nonconvex_hull = FALSE)

## Add covariates to the BAUs
d_cutoff <- 30.48  
d_BAU <- distR(coordinates(BAUs), GZ_df)
BAUs$x1 <- as.numeric(d_BAU < d_cutoff)
BAUs$x2 <- d_BAU * BAUs$x1
BAUs$x3 <- as.numeric(d_BAU >= d_cutoff)
BAUs$x4 <- d_BAU * (BAUs$x3)


# ---- FRK analysis ---- 

cat("Starting FRK analysis...\n")

## In FRK, covariates are supplied at the BAU level: Remove covariates from the 
## data as the BAUs and the data cannot have common fields
Am_data@data[, c("x1", "x2", "x3", "x4")] <- NULL

# coordinates(Am_data) = ~ Easting + Northing

BAUs$fs <- 1     # scalar matrix for fine scale variation
Am_data$std <- 1 # set measurement error to small value to replicate lognormal kriging

suppressWarnings( # Suppress warning about log-link being inappropriate for Gaussian data 
  M <- FRK(f = Am ~ -1 + x1 + x2 + x3 + x4, data = list(Am_data),
           response = "gaussian", link = "log", BAUs = BAUs, est_error = FALSE, nres = 2,
           # manually set these arguments to reduce console output:
           K_type = "precision", method = "TMB")
)


## Continuous predictions
pred      <- predict(M, type = c("link", "mean"))
plot_list <- plot(M, pred$newdata)

formatter <- function(x){ 
  x/10^5 
}
x_scale <- scale_x_continuous(breaks = 10^5 *c(2.197, 2.199, 2.201), labels = formatter)  
y_scale <- scale_y_continuous(breaks = 10^5 * c(2.852, 2.854, 2.856), labels = formatter)
plot_list <- lapply(
  plot_list, 
  function(gg) gg + lab1 + lab2 + x_scale + y_scale)

figure <- ggpubr::ggarrange(
  plot_list$p_Y + labs(fill = expression(bold(Y)~"pred.")) +   
    scale_fill_gradientn(colours = nasa_palette, labels = scales::scientific),
  plot_list$RMSPE_Y + labs(fill = expression(RMSPE(bold(Y)))),
  plot_list$p_mu + labs(fill = expression(bold(mu)~"pred.")) +    
    scale_fill_gradientn(colours = nasa_palette, labels = scales::scientific),
  plot_list$RMSPE_mu + labs(fill = expression(RMSPE(bold(mu)))), 
  align = "hv", nrow = 2, ncol =2) 

ggsave( 
  figure,
  filename = "4_2_Am_FRK_cts.png", device = "png", width = 9.6, height = 5.7,
  path = "Figures/"
)


## Predict over the blocks 
RNGversion("3.6.0"); set.seed(1)
pred <- predict(M, newdata = blocks, percentiles = NULL)

FRK_results <- pred$newdata@data %>% 
  mutate(area_sqrt = sapply(blocks@polygons, function(x) sqrt(x@area)), 
         Scheme = as.numeric(blocks@data$Scheme), 
         Framework = "FRK") %>%
  dplyr::select(c("p_mu", "RMSPE_mu", "area_sqrt", "Scheme", "Framework")) %>% 
  melt(id.vars = c("area_sqrt", "Scheme", "Framework"))

cat("FRK analysis complete.\n")


# ---- Only FRK ----

df <- FRK_results %>% 
  dplyr::select(c("variable", "value", "area_sqrt", "Scheme")) %>% 
  # alter the labels to change the facet_wrap titles:
  mutate(
    variable = factor(
    variable, 
    labels = c("'Block prediction'", "'RMSPE from block prediction'")
    )) %>% 
  mutate(
    Scheme = as.character(Scheme)
  )


figure <- ggplot(data = df, aes(x = area_sqrt)) +
  geom_line(aes(y = value), size = 1) +
  facet_wrap(~variable, scales = "free", labeller = label_parsed) + 
  labs(x = "block size (m)", y = "") +
  theme_bw() + 
  scale_y_continuous(labels = scales::scientific) + 
  theme(text = element_text(size = 20), 
        strip.text = element_text(size = 20))

figure <- figure + 
  theme(
    strip.background = element_blank()
  )

ggsave( 
  figure,
  filename = "4_2_Am_FRK.png", device = "png", width = 13.6, height = 4.5,
  path = "Figures/"
)

