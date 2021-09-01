library("FRK") 
library("dplyr")
library("ggplot2")
library("ggmap")
library("spacetime")
library("maptools") # readShapePoly()
library("reshape2")
library("sp")
library("stringr")
library("htmltab")
library("ggpubr")

# ---- Preprocessing ----

source("./src/Chicago_prep.R")


# ---- BAUs ----

## Set up the space-time BAUs. 
## For the spatial BAUs, we use the Chicago community areas.
ST_BAUs <- auto_BAUs(manifold = STplane(),
                     data = chicago_crimes_fit,
                     spatial_BAUs = community_areas,
                     tunit = "years")
ST_BAUs$fs <- 1 # scalar fine-scale variance matrix

## Create population covariate
ST_BAUs$population <- community_areas$population

## Add covariates to BAUs that can be used to make a piecewise-linear temporal trend
## As the inclusion of the trend did not significantly affect the results, it 
## was omitted from the paper to simplify the exposition.

# year <- ST_BAUs@data$t + 2000
# ST_BAUs$x1 <- as.numeric(year < 2014)
# ST_BAUs$x2 <- year * ST_BAUs$x1
# ST_BAUs$x3 <- as.numeric(year >= 2014)
# ST_BAUs$x4 <- year * ST_BAUs$x3

# ---- Model fitting ----

basis <- auto_basis(STplane(), chicago_crimes_fit, tunit = "years")

M <- FRK(f = number_of_crimes ~ log(population),   
         data = chicago_crimes_fit, basis = basis, BAUs = ST_BAUs,         
         response = "poisson", link = "log", 
         sum_variables = "number_of_crimes", fs_by_spatial_BAU = TRUE) 

print(object.size(M), units = "Mb")
Chicago_SRE_object <- M
saveRDS(Chicago_SRE_object, file = "./intermediates/Chicago_SRE_object.rds")

# ---- Simplified, high-level version that I show in presentations ---- 

# M <- FRK(f = number_of_crimes ~ 1,
#          response = "poisson",
#          link     = "log",
#          data     = chicago_crimes,
#          sum_variables = "number_of_crimes")

# ---- Prediction ----

RNGversion("3.6.0")
set.seed(1)
pred <- predict(M, type = "response",  
                percentiles = c(5, 95, 10, 90, 15, 85, 20, 80, 25, 75))  

## Shift t by 2000, so we can refer to it by the year rather than temporal index
## (NB: it just so happens that the first year in this dataset is 2001)
pred$newdata@data$t <- pred$newdata$t + 2000 

## STFDF with predictions and prediction uncertainty, which we will add to
ST_pred <- pred$newdata


# ---- Compute binned validation data ----

## Recall that in the model fitting stage, we left out the data from the years 
## 2010 and 2019. To perform model validation on our BAU level predictions in 
## these years, we first need to obtained the binned data at the BAUs for 
## these years. The simplest approach is to simply map the entire dataset to 
## the BAUs, and then subset the data corresponding to validation years.

## Map the prediction and forecast crime data to BAUs. 
binned_data <- FRK:::map_data_to_BAUs(chicago_crimes, sp_pols = ST_pred, 
                                      sum_variables = "number_of_crimes")

## Extract the data and arrange by the BAU id
df_val <- binned_data@data %>% arrange(n)

## Assign the observed count to the prediction data frame for plotting
ST_pred@data$number_of_crimes <- df_val$number_of_crimes


# ---- Prediction and forecasting years ----

## Wrap the plotting code in a rudimentary function for convenience
## NB: subset_time are the years we wish to visualise
plot_predictions <- function(subset_time) {
  
  ## Plot the predictions and uncertainty
  plots <- plot(M, pred$newdata,
                map_layer = chicago_map, subset_time = subset_time, 
                colour = "black", size = 0.3, alpha = 0.85)
  
  ## plot the validation data. 
  ## NB: plot() does plot the data used in model fitting, but since 2010 and 
  ## 2019 are validation years (i.e., those years were omitted during model 
  ## fitting), we don't have any data to plot! 
  ## Construct the plot manually:
  plots$number_of_crimes <- plot_spatial_or_ST(
    ST_pred, all.vars(M@f)[1],
    map_layer = chicago_map, subset_time = subset_time, 
    colour = "black", size = 0.3, alpha = 0.85
  )[[1]]
  
  ## Change layout of each quantity to a single column, and edit axis
  plots <- lapply(plots, function(gg) gg + 
                    facet_wrap(~t,  ncol = 1) + 
                    xlab("lon (deg)") + ylab("lat (deg)"))
  
  ## Ensure predictions and observed are on same scale, for easy comparisons
  count_lims <- ST_pred@data %>% 
    subset(t %in% c(2010, 2019)) %>%
    select(c("number_of_crimes", "p_Z")) %>%
    c() %>% range()
  
  ## Edit titles and legends
  plots$number_of_crimes <- plots$number_of_crimes + 
    scale_fill_distiller(palette = "Spectral",  breaks = c(1000, 3000, 5000), lim = count_lims) + 
    labs(title = "Observed (withheld)\ncrime count", fill = "") + 
    theme(axis.title.x = element_blank())
  
  plots$p_Z <- plots$p_Z + 
    scale_fill_distiller(palette = "Spectral", breaks = c(1000, 3000, 5000), lim = count_lims) + 
    labs(title = "Predicted \ncrime count", fill = "") + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank())
  
  plots$interval90_Z <- plots$interval90_Z + 
    scale_fill_distiller(palette = "BrBG", n.breaks = 3) + 
    labs(title = "90% prediction-\ninterval width", fill = "") + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.title.x = element_blank(), 
          axis.ticks.y = element_blank())
  
  ## Edit visual aspects of the plot, such a font size and legend details
  plots <- lapply(plots, function(gg) gg + 
                    theme(axis.text = element_text(size = 12),
                          axis.title = element_text(size = 14), 
                          legend.text = element_text(size = 12), 
                          plot.title = element_text(hjust = 0.5, size = 14), 
                          legend.key.width = unit(0.9, "cm"), 
                          strip.text = element_text(size = 12), 
                          plot.margin = unit(c(0,0,0,0), "lines")) +
                    scale_x_continuous(breaks = c(-87.6, -87.8), expand = c(0, 0)) + 
                    scale_y_continuous(n.breaks = 3, expand = c(0, 0)))
  
  return(plots)
}

## First plot both the pred and forecast years
subset_time <- c(10, 19)
plots <- plot_predictions(subset_time = subset_time)
figure <- ggarrange(plots$number_of_crimes, plots$p_Z, plots$interval90_Z, 
                    align = "hv", nrow = 1, legend = "top")

ggsave(figure, filename = "Chicago_data_pred_uncertainty.png", 
       device = "png", width = 10, height = 8.5, path = "./img/")

## Also save a version with only the forecast years (for presentations)
subset_time <- 19
plots <- plot_predictions(subset_time = subset_time) 
figure <- ggarrange(plots$number_of_crimes, plots$p_Z, plots$interval90_Z, 
                    align = "hv", nrow = 1, legend = "top")
ggsave(figure, filename = "Chicago_data_pred_uncertainty_2019_only.png", 
       device = "png", width = 10, height = 5.5, path = "./img/")


# ---- Time series plot ----

## Now a time-series of three community areas of interest, showing at each 
## time-point the prediction, true observed value, and the prediction interval. 

RNGversion("3.6.0")
set.seed(1996)
ns <- length(ST_BAUs@sp)
spBAU_id <- sample(1:ns, 3, replace = FALSE)
cat("Focused community areas: ", 
    paste(community_areas$community[spBAU_id], collapse=", "), "\n")

nt <- length(ST_BAUs@time)
increment  <- ns * 0:(nt - 1) # Terms to add to each rmidx
ST_BAU_id  <- sort(c(outer(spBAU_id, increment, "+")))

ST_BAU_community_name <- as.character(community_areas@data$community[rep(spBAU_id, times = nt)])

time_series_df <- ST_pred@data[ST_BAU_id, ] %>% 
  mutate(community = str_to_title(ST_BAU_community_name))

## Create a long form version for ggplot (so we can split prediction and observed by colour)
time_series_df_long <- time_series_df %>% 
  select(community, t, number_of_crimes, p_Z) %>% 
  melt(id = c("community", "t"))

time_series_plot <- ggplot() + 
  geom_vline(xintercept = c(2010, 2019), colour = "grey", alpha = 0.5, size = 1.5) +
  geom_errorbar(data = time_series_df, 
                aes(x = t, y = p_Z, ymin = Z_percentile_5, ymax = Z_percentile_95),
                width = 0.25, alpha = 0.5) +
  geom_point(data = time_series_df_long, aes(x = t, y = value, colour = variable), 
             size = 0.9 * (time_series_df_long$variable == "p_Z") + 1.2 * (time_series_df_long$variable == "number_of_crimes")) +
  facet_wrap(~community, nrow = 3, ncol = 1, scales = "free_y") +
  labs(colour = "", x = "Year", y = "Number of crimes") +
  scale_colour_discrete(labels = c("Observed count", "Predicted count"))+
  theme_bw() + 
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 13), 
        legend.text = element_text(size = 13), 
        strip.text = element_text(size = 12), 
        legend.position = "top") + 
  scale_x_continuous(breaks = c(2001, 2005, 2010, 2015, 2019)) + 
  scale_y_continuous(n.breaks = 4)

ggsave(time_series_plot,
  path = "./img/", filename = "Chicago_focused_CAs_time_series.png", 
  device = "png", width = 9, height = 6)

## Select the MC samples corresponding to our chosen BAUs
n_MC <- ncol(pred$MC$Z_samples)
MC_df <- data.frame(community = rep(time_series_df$community, each = n_MC),
                    t = ST_pred$t %>%
                      unique() %>%
                      rep(each = length(unique(time_series_df$community)) * n_MC) %>%
                      as.factor(),
                    samples = as.vector(t(pred$MC$Z_samples[ST_BAU_id, ])))

MC_df_val <- MC_df %>% subset(t %in% c(2010, 2019))
time_series_df_long_val <- time_series_df_long  %>% subset(t %in% c(2010, 2019))
MC_df_val$t <- MC_df_val$t %>% droplevels()

predictive_distribution_plots <- ggplot() + 
  geom_density(data = MC_df_val, 
               aes(x = samples)) +
  geom_segment(data = time_series_df_long_val, 
               aes(x = value, xend = value, y = 0, yend = Inf, colour = variable), 
               size = 0.5 * (time_series_df_long_val$variable == "p_Z") + 1 * (time_series_df_long_val$variable == "number_of_crimes")) +
  facet_grid(t ~ community, scales = "free_x") + 
  labs(colour = "", x = "Number of crimes") +
  scale_colour_discrete(labels = c("Observed count", "Predicted count")) + 
  theme_bw() + 
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 13), 
        legend.text = element_text(size = 13), 
        strip.text = element_text(size = 12), 
        legend.position = "top") 

ggsave(predictive_distribution_plots,
  path = "./img/", filename = "Chicago_focused_CAs_predictive_distributions.png", 
  device = "png", width = 9, height = 6)

# ---- Coverage ----

## Note that we will use the response prediction intervals, as the prediction 
## intervals of the mean mu will NEVER contain zero, so it would be impossible 
## for BAUs with zero counts to have their data in the prediction interval. 

## Empirical coverage, and the MAPE (how close was the predicted crime to 
## the true crime in each year (on average)?
Chicago_coverage_and_MAPE <- ST_pred@data %>% 
  subset(t %in% c(2010, 2019)) %>% 
  mutate(
    data_in_pred_interval_90 = Z_percentile_5 <= number_of_crimes & number_of_crimes <= Z_percentile_95,
    data_in_pred_interval_80 = Z_percentile_10 <= number_of_crimes & number_of_crimes <= Z_percentile_90,
    data_in_pred_interval_70 = Z_percentile_15 <= number_of_crimes & number_of_crimes <= Z_percentile_85,
    data_in_pred_interval_60 = Z_percentile_20 <= number_of_crimes & number_of_crimes <= Z_percentile_80,
    absolute_percentage_error = abs(p_Z - number_of_crimes)/number_of_crimes
  ) %>%
  group_by(t) %>%
  summarise(
    ## Also multiply by 100 so that reported values are a percentage
    coverage_90 = mean(data_in_pred_interval_90) * 100,
    coverage_80 = mean(data_in_pred_interval_80) * 100,
    coverage_70 = mean(data_in_pred_interval_70) * 100,
    coverage_60 = mean(data_in_pred_interval_60) * 100, 
    MAPE = mean(absolute_percentage_error) * 100
  ) %>% 
  as.data.frame()


Chicago_coverage_and_MAPE$average_coverage_difference <- 
  rowMeans(Chicago_coverage_and_MAPE[, paste0("coverage_", seq(90, 60, -10))] - 
             matrix(rep(seq(90, 60, -10), 2), nrow = 2, byrow = TRUE))

write.csv(Chicago_coverage_and_MAPE, 
          "./results/Chicago_coverage_and_MAPE.csv", 
          row.names = FALSE)

