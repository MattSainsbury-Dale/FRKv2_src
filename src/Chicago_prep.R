# ---- Chicago ggmap ----

load("./data/chicago_map.RData")
load("./data/chicago_bbox.RData")

## Create map layer to place under all plots
chicago_map <- ggmap(chicago)

## Define x (longitude) and y (latitude) limits here. This will restrict the plots
## that use chicago_map as a base layer. 
x_lims <- xlim(c(chicago_bbox["left"], chicago_bbox["right"])) 
y_lims <- xlim(c(chicago_bbox["bottom"], chicago_bbox["top"])) 

## Add axis labels/limits and theme choices to chicago_map to reduce code repetition
chicago_map <- chicago_map + 
  x_lims + y_lims +
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
tab <- htmltab("https://en.wikipedia.org/w/index.php?title=Community_areas_in_Chicago&oldid=1012989271",1)
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
  dplyr::summarise(number_of_crimes = n()) %>%
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

