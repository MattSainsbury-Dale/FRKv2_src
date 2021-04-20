library("FRK")
library("sp")
library("dplyr")
library("reshape2")

Am_data <- read.csv("./intermediates/Am_data.csv")
BAUs <- readRDS("./intermediates/Am_BAUs.rds")
blocks <- readRDS("./intermediates/Am_blocks.rds")


# ---- FRK ----

## In FRK, covariates are supplied at the BAU level: Remove covariates from the 
## data as the BAUs and the data cannot have common fields
Am_data[, c("x1", "x2", "x3", "x4")] <- NULL

coordinates(Am_data) = ~ Easting + Northing

BAUs$fs <- 1 # scalar matrix for fine scale variation
Am_data$std <- 1 # set measurement error to small value to replicate lognormal kriging

# M <- FRK(f = Am ~ x1 + x2 + x3, data = list(Am_data),
M <- FRK(f = Am ~ -1 + x1 + x2 + x3 + x4, data = list(Am_data),
         response = "gaussian", link = "log",
         BAUs = BAUs, est_error = FALSE)


# # ---- Predict and plot over the BAUs ----
# 
# RNGversion("3.6.0")
# set.seed(1)
# pred <- predict(M, type = c("link", "mean"), percentiles = NULL)
# plotlist <- plot(M, pred)
# 
# ggsave( 
#   ggarrange(plotlist = plotlist),
#   filename = "Americium_BAU_predictions.png", device = "png", width = 9.5, height = 4.5,
#   path = "./img/"
# )

# ---- Predict over the blocks ----

## Now we predict over arbitrary polygons.

## Predict using FRK
RNGversion("3.6.0")
set.seed(1)
pred <- predict(M, newdata = blocks, percentiles = NULL)

FRK_results <- pred$newdata@data
FRK_results$area_sqrt <- sapply(blocks@polygons, function(x) sqrt(x@area))
FRK_results$Scheme <- as.numeric(blocks@data$Scheme)
FRK_results$Framework <- "FRK" 

FRK_results <- FRK_results %>%
  dplyr::select(c("p_mu", "RMSPE_mu", "area_sqrt", "Scheme", "Framework")) %>% 
  melt(id.vars = c("area_sqrt", "Scheme", "Framework"))

write.csv(FRK_results, 
          file = "./intermediates/Am_FRK.csv", 
          row.names = FALSE)



# ## Plotting over blocks is a bit tricky... just do it manually rather 
# ## than using plot()
# block_pred_df <- FRK:::.SpatialPolygonsDataFrame_to_df(pred$newdata)
# 
# ## Change level order to reverse the order that the blocks are plotted.
# n_schemes <- 2; n_block <- 5
# block_pred_df$block <- factor(block_pred_df$id,
#                               levels = paste0("block", (n_schemes * n_block):1))
# p_block <- p_basic +
#   geom_point(size = 0.3) +
#   geom_polygon(data = block_pred_df, aes(fill = p_mu, group = block),
#                alpha = 1, colour = "black") +
#   scale_fill_distiller(
#     palette = "Spectral",
#     name = expression(paste(widehat(p)[mu]["|"][bold(Z)], " \U2261 ", "E(", mu *"(\U00B7)", " | ", bold(Z), ", ", bold("\u03B8"), ") ")), 
#     limits = range(Am_data$Am), 
#     labels = scales::scientific
#     )
# 
# p_block_RMSPE <- p_basic +
#   geom_point(size = 0.3) +
#   geom_polygon(data = block_pred_df, aes(fill = RMSPE_mu, group = block),
#                alpha = 1, colour = "black") +
#   scale_fill_distiller(palette = "BrBG",
#                        name = expression(RMSPE(widehat(p)[mu]["|"][bold(Z)], mu)),
#                        labels = scales::scientific)
# 
# ggsave(
#   ggpubr::ggarrange(p_block, p_block_RMSPE, nrow = 1, align = "hv"),
#   filename = "Americium_block_predictions.png", device = "png", width = 9.5, height = 4.5,
#   path = "./img/"
# )



