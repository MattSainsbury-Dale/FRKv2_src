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


defaultW <- getOption("warn") 
options(warn = -1)


# ---- Chicago ggmap ----

load("./data/chicago_map.RData")

## Create map layer to place under all plots
chicago_map <- ggmap(chicago)

## Define x (longitude) and y (latitude) limits here, so all plots are consistent
x_lims <- xlim(c(-87.95,-87.49)) 
y_lims <- ylim(c(41.62, 42.04))

## Add axis labels/limits and theme choices to chicago_map to reduce code repetition
chicago_map <- chicago_map + x_lims + y_lims +
  xlab("lon (deg)") + ylab("lat (deg)") + 
  theme_bw() + coord_fixed()


# ---- Chicago crime dataset ----

## Full list of crimes in Chicago between 2000 and 2019
load("./data/chicago_crime_df.Rda")

## There are some invalid locations. Just remove these observations:
df <- subset(df, location != "")

## After plotting, I found some locations which are way outside the city.
## We could remove anything south of latitude 41 and west of longitude -89: 
## df <- subset(df, latitude > 41 & longitude > -89)
## However, it's not necessary because FRK will just exclude observations that 
## did not fall in the community areas. Removal would be necessary if we were 
## using auto_BAUs() for the spatial BAUs. 

## We focus on violent, non-sexual crimes. These crimes include "ASSAULT", 
## "BATTERY" and "HOMICIDE". The frequency of Assult or Battery is much higher
## than that of homicide (as one would expect). Furthermore, homicdie is a much 
## more serious crime than assault or battery. Hence, we will omit homicide 
## and focus only on assault and battery.
# sum(df$primary_type == "HOMICIDE")
df <- subset(df, primary_type %in% c("ASSAULT", "BATTERY"))
df <- droplevels(df) # Drop unused levels



# ---- Chicago community areas ----

## Load Shapefile of community areas, and name the coordinates
community_areas <- readShapePoly("./data/Chicago_shapefiles/chicago_community_areas.shp")
coordnames(community_areas) <- c("longitude", "latitude")

## The polygon IDs are labelled from 0 to 76, and are not even associated 
## with the community number:
# sapply(community_areas@polygons, slot, "ID") 
# community_areas@data$area_num_1
## This could cause confusion later, so lets rename the polygons now.
new_IDs <- as.character(community_areas@data$area_num_1)
for (i in 1:length(community_areas@polygons)){
  slot(community_areas@polygons[[i]], "ID") = new_IDs[i]
}

## Sanity check:
# all(sapply(community_areas@polygons, slot, "ID")  == as.character(community_areas@data$area_num_1))

## Remove the community area with the airport (O'Hare):
community_areas <- subset(community_areas, area_num_1 != 76)

## We removed one of the areas: remove it from the levels
community_areas@data <- droplevels(community_areas@data)

## Reset the row names
rownames(community_areas@data) <- NULL


# ---- Population covariates (community area level) ----

## Read the HTML table on the community areas wikipedia page 
tab <- htmltab("https://en.wikipedia.org/wiki/Community_areas_in_Chicago",1)
rownames(tab) <- NULL

## The last row corresponds to column totals, so remove it
tab <- tab[-nrow(tab), ]

## The table is in ascending order in terms of community area (i.e., 
## it is ordered from 1 to 77). Use the area number in the shapefile object
## to map the population. Note that O'Hare will simply not be selected.
tmp <- tab[, "2017 population"]
tmp <- gsub(",", "", tmp)
tmp <- as.numeric(tmp)
idx <- community_areas@data$area_num_1 %>% as.character %>% as.numeric
## (note that we needed to use as.character before as.numeric. This is 
## because factors are stored internally as integers with a table to 
## give the factor level labels. Just using as.numeric will only give 
## the internal integer codes.)
community_areas$population <- tmp[idx]

g_population <- plot_spatial_or_ST(community_areas, "population", map_layer =  chicago_map,  
           colour = "black", size = 0.3)[[1]] + 
  geom_text(data = cbind(data.frame(community_areas), coordinates(community_areas)), 
            aes(label = area_num_1), size = 2)

# ---- Spatio-temporal dataframe (one-year time periods) ----

## Create spatio-temporal dataframe.
## Given that 2020 is incomplete, we will omit it from this study. 
ST_df <- df %>% 
  subset(year != 2020) %>%
  group_by(year, longitude, latitude) %>%
  summarise(number_of_crimes = n()) %>%
  as.data.frame()

## Visualization: Number of crimes in each year:
g_temporal_trend <- ggplot(data = ST_df %>%
         group_by(year) %>%
         summarise(total_crimes = sum(number_of_crimes)), 
       aes(x = year, y = total_crimes)) +
  geom_point() +
  geom_smooth(colour = "red", method = 'lm', se = F) + 
  geom_smooth(aes(group = year < 2014), method = 'lm', se = F) + 
  labs(y = "total crimes") + 
  theme_bw()
## A piecewise temporal trend seems appropriate

## Create a Date field
ST_df$time <- as.Date(paste(ST_df$year, 06,01,sep="-"))

## Construct an STIDF object
chicago_crimes <- stConstruct(x = ST_df,                               
                              space = c("longitude", "latitude"), 
                              time = "time",                      
                              interval = TRUE)      # time reflects an interval

## We leave the years 2010 and 2019 for model validation.
chicago_crimes_fit <- subset(chicago_crimes, 
                             !(chicago_crimes@data$year %in% c(2010, 2019, 2020)))


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

## Add covariates to BAUs that can be used to make a piecewise linear trend
year <- ST_BAUs@data$t + 2000
ST_BAUs$x1 <- as.numeric(year < 2014)
ST_BAUs$x2 <- year * ST_BAUs$x1
ST_BAUs$x3 <- as.numeric(year >= 2014)
ST_BAUs$x4 <- year * ST_BAUs$x3

# ---- Basis functions ----

basis <- auto_basis(STplane(), chicago_crimes_fit, tunit = "years", nres = 3)

# ---- Model fitting ----

M <- FRK(f = number_of_crimes ~ -1 + sqrt(population) + x1 + x2 + x3 + x4,   
         data = list(chicago_crimes_fit), basis = basis, BAUs = ST_BAUs,         
         response = "poisson", link = "square-root", 
         sum_variables = "number_of_crimes", fs_by_spatial_BAU = TRUE) 

print(object.size(M), units = "Mb")
Chicago_SRE_object <- M
saveRDS(Chicago_SRE_object, file = "./intermediates/Chicago_SRE_object.rds")

# ---- Prediction ----

RNGversion("3.6.0")
set.seed(1)
system.time(
  pred <- predict(M,type = "response", 
                  percentiles = c(5, 95, 10, 90, 15,85, 20, 80,25,75))  
)


## Shift t by 2000, so we can refer to it by the year rather than temporal index
## (NB: it just so happens that the first year in this dataset is 2001)
pred$newdata@data$t <- pred$newdata$t + 2000 

## STFDF with predictions and prediction uncertainty, which we will add to
ST_pred <- pred$newdata


# ---- Compute binned validation data ----

## Recall that in the model fitting stage, we left out the data from the years 
## 2010 and 2019. To perform model validation on our BAU level predictions in 2010
## and 2019, we first need to obtained the binned data at the BAUs for these years.
## The simplest approach is to simply map the entire dataset to the BAUs, 
## and then subset the data corresponding to validation years.

## Map the prediction and forecast crime data to BAUs. 
binned_data <- FRK:::map_data_to_BAUs(chicago_crimes, sp_pols = ST_pred, 
                                      sum_variables = "number_of_crimes")

## Extract the data and arrange by the BAU id
df_val <- binned_data@data %>% arrange(n)

## Sanity checks: 
## ensure the the observed counts are in the same order as the predictions.
all(ST_pred$latitude == df_val$latitude)
all(ST_pred$longitude == df_val$longitude)
all(ST_pred$year == df_val$year)
all(ST_pred$n == df_val$n)

## Assign the observed count to the prediction data frame for plotting
ST_pred@data$number_of_crimes <- df_val$number_of_crimes


# ---- Prediction and forecasting years ----

subset_time <- c(2010, 2019) ## Years we wish to analyse

## NB: plotting zdf works provided the column name to plot is all.vars(M@f)[1]
plots <- plot(M, pred, zdf = ST_pred, 
              map_layer = chicago_map, subset_time = subset_time, 
              colour = "black", size = 0.3, alpha = 0.85)

## Change layout of each quantity to a single column, and edit axis
plots <- lapply(plots, function(gg) gg + 
                  facet_wrap(~t,  ncol = 1) + 
                  xlab("lon (deg)") + ylab("lat (deg)"))

## Ensure predictions and observed are on same scale, for easy comparisons
count_lims <- ST_pred@data %>% 
  subset(t %in% subset_time) %>%
  select(c("number_of_crimes", "p_Z")) %>%
  c() %>%range()

suppressMessages(ggsave( 
  ggarrange(
    plots$number_of_crimes + scale_fill_distiller(palette = "Spectral", n.breaks = 4, name = "Observed \ncrime count", lim = count_lims), 
    plots$p_Z  +  scale_fill_distiller(palette = "Spectral", n.breaks = 4, name = "Predicted \ncrime count", lim = count_lims), 
    plots$interval90_Z + scale_fill_distiller(palette = "BrBG", n.breaks = 4, name = "90% prediction \ninterval width"), 
    align = "hv", nrow = 1, legend = "top"),
  filename = "Chicago_data_pred_uncertainty.png", device = "png", width = 11, height = 7,
  path = "./img/"
))

# ---- Time series plot ----


## Now a time-series of three community areas of interest, showing at each time-point
## the prediction, true observed value, and the prediction interval. 

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

suppressMessages(ggsave( 
  ggplot() + 
    geom_vline(xintercept = c(2010, 2019), colour = "grey", alpha = 0.5, size = 1.5) +
    geom_errorbar(data = time_series_df, 
                  aes(x = t, y = p_Z, ymin = Z_percentile_5, ymax = Z_percentile_95),
                  width = 0.25, alpha = 0.5) +
    geom_point(data = time_series_df_long, aes(x = t, y = value, colour = variable), 
               size = 0.9 * (time_series_df_long$variable == "p_Z") + 1.2 * (time_series_df_long$variable == "number_of_crimes")) +
    facet_wrap(~community, nrow = 3, ncol = 1, scales = "free_y") +
    labs(colour = "", x = "Year", y = "Number of crimes") +
    scale_colour_discrete(labels = c("Observed count", "Predicted count"))+
    theme_bw()
  ,
  filename = "Chicago_focused_CAs_time_series.png", device = "png", width = 9, height = 6,
  path = "./img/"
))

# ## Using geom_line() and geom_ribbon() for the shaded area (not appropriate for highly discretised time)
# ggplot() + 
#   geom_line(data = time_series_df_long, aes(x = t, y = value, colour = variable)) + 
#   geom_ribbon(data = time_series_df, 
#               aes(x = t, ymin = Z_percentile_5, ymax = Z_percentile_95), 
#               fill= "gray", alpha = 0.5) +
#   facet_wrap(~community, nrow = 3, ncol = 1, scales = "free_y") + 
#   labs(colour = "", x = "Year", y = "Number of crimes") +
#   scale_colour_discrete(labels = c("True data", "Prediction")) + 
#   theme_bw()

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

suppressMessages(ggsave( 
  ggplot() + 
    geom_density(data = MC_df_val, 
                 aes(x = samples)) +
    geom_segment(data = time_series_df_long_val, 
                 aes(x = value, xend = value, y = 0, yend = Inf, colour = variable), 
                 size = 0.5 * (time_series_df_long_val$variable == "p_Z") + 1 * (time_series_df_long_val$variable == "number_of_crimes")) +
    facet_grid(t ~ community, scales = "free_x") + 
    labs(colour = "", x = "Number of crimes") +
    scale_colour_discrete(labels = c("Observed count", "Predicted count")) + 
    theme_bw()
  ,
  filename = "Chicago_focused_CAs_predictive_distributions.png", device = "png", width = 9, height = 6,
  path = "./img/"
))

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

options(warn = defaultW)
