# ---- Census data pre-processing ----

## Load the data
census_SA1_df <- read.csv("./data/Sydney_SA1_data.csv")
census_SA2_df <- read.csv("./data/Sydney_SA2_data.csv")

## Poverty lines as defined in Appendix "Sydney poverty lines"
poverty_lines <- c("Couple_family_with_no_children" = round_any(594.6, 200), 
                   "Couple_family_with_children" = round_any(835.3, 200), 
                   "One_parent_family" = round_any(691.04, 200))

## Select the income brackets common to all family groups as defined above
low_income <- c("Negative_Nil_income_", "X1_199_", "X200_299_", 
                "X300_399_", "X400_599_")

## Select the family names we are interested in.
family_names <- names(poverty_lines)

## The following function performs data pre-processing on the SA1 or SA2 data frames.
## At each SA, it first determines the number of families below the poverty line 
## within each family type, then computes the total number of families in poverty
## (i.e., the number of "successes", Z) and the number of families in total 
## (i.e., the size parameter, k). For the sake of data visualisation, we also 
## compute the empirical proportion of families in poverty.
censusDataPreprocess <- function(df, SA) {
  
  df <- df %>% 
    mutate(
      poverty_count_Couple_family_with_no_children = rowSums(.[, paste0(low_income, family_names[1])]),
      poverty_count_Couple_family_with_children = rowSums(.[, paste0(low_income, family_names[2])]), 
      poverty_count_One_parent_family = rowSums(.[, paste0(c(low_income, "X600_799_"), family_names[3])])
      
    ) %>% 
    mutate(
      total_poverty_count = poverty_count_Couple_family_with_no_children + 
        poverty_count_Couple_family_with_children + 
        poverty_count_One_parent_family
    ) %>% 
    mutate(
      Total_families_of_interest = 
        Total_Couple_family_with_no_children + 
        Total_Couple_family_with_children + 
        Total_One_parent_family
    )  %>% 
    mutate(
      Proportion_poverty = total_poverty_count / Total_families_of_interest
    ) %>% 
    dplyr::select(region_id, total_poverty_count, Total_families_of_interest, Proportion_poverty)
  
  return(df)
}

## Compute the total number of families, number of families in poverty, 
## and the proportion of families in poverty at SA1 and SA2 level. 
census_SA1_df <- censusDataPreprocess(census_SA1_df, "SA1")
census_SA2_df <- censusDataPreprocess(census_SA2_df, "SA2")


# ---- SA Shapefiles and pre-processing ----

## Since SA3s are aggregations of SA2s, which in turn are
## aggregations of SA1s, we define the domain of interest first in terms of 
## SA3s, and then subset the SA2s and SA1s accordingly. This is to ensure that 
## the prediction maps look consistent, in particular, the SA3 predictions don't 
## have a big chunk of SA1s missing.  This will ensure that the SA1s and SA2s 
## used for model fitting aggregate to whole SA3s. 

## Deal with the SA3s first so that we can subset the SA1s based on the SA3s. 
suppressWarnings(
  SA3 <- readShapePoly("./data/Sydney_shapefiles/SA3/SA3_2011_AUST.shp", delete_null_obj = TRUE)
)
coordnames(SA3) <- c("lon", "lat")

## Now focus on a subset of SA3s in NSW (around Sydney)
cc          <- coordinates(SA3)
idx         <- cc[, "lon"] > 150.77 & cc[, "lon"] < 151.5 & cc[, "lat"] > -34.1 & cc[, "lat"] < -33.73
SA3_NSW_sub <- SA3[idx,]
# number of SA3s used for fitting
cat("Number of SA3s:", length(SA3_NSW_sub), "\n")

## Deal with the SA2s first so that we can subset the SA1s based on the SA2s 
SA2 <- readShapePoly("./data/Sydney_shapefiles/SA2/SA2_2011_AUST.shp", delete_null_obj = TRUE)
coordnames(SA2) <- c("lon", "lat")

## Add the SA2 census data 
census_SA2_df$SA2_MAIN11 <- census_SA2_df$region_id
SA2@data <- merge(SA2@data, census_SA2_df, by = "SA2_MAIN11")

## Retain only those SA2s which are in the remaining SA3s
SA2_NSW_sub <- subset(SA2, SA3_CODE11 %in% SA3_NSW_sub$SA3_CODE11) 

## Remove any SA2s which had 0 total families of interest
SA2_NSW_sub_withk0 <- SA2_NSW_sub  # keep a backup with all of the SA2s
idx <- which(SA2_NSW_sub@data$Total_families_of_interest == 0)
cat("SA2s with no families of interest: ", 
    paste(SA2_NSW_sub@data$SA2_NAME11[idx], collapse=", "), "\n")
SA2_NSW_sub <- subset(SA2_NSW_sub, Total_families_of_interest > 0)
cat("Number of SA2s with at least 1 family of interest:", length(SA2_NSW_sub), "\n")

## Now deal with the SA1s
SA1 <- readShapePoly("./data/Sydney_shapefiles/SA1/SA1_2011_AUST.shp", delete_null_obj=TRUE)
coordnames(SA1) <- c("lon", "lat")

## Retain only those SA1s which are in the remaining SA2s.
SA1_NSW_sub <- subset(SA1, SA2_MAIN11 %in% SA2_NSW_sub$SA2_MAIN11) 
cat("Number of SA1s:", length(SA1_NSW_sub), "\n")

## Add the SA1 census data
census_SA1_df$SA1_7DIG11 <- census_SA1_df$region_id
SA1_NSW_sub@data <- merge(SA1_NSW_sub@data, census_SA1_df,
                          by = "SA1_7DIG11")

## For code simplicty, remove _NSW_sub henceforth
SA1s <- SA1_NSW_sub
SA2s <- SA2_NSW_sub


construct_training_data <- function(fitting = "mixed", SA1s, SA2s) {
  if (fitting == "SA2s") {
    return(SA2s)
  } else if (fitting == "SA1s") {
    return(SA1s)
  } else if (fitting == "mixed") {
    
    ## We include some SA1s in the training set 
    SA1s <- SA1s
    zdf  <- SA2s 
    
    ## Select some data supports to be replaced
    RNGversion("3.6.0")
    set.seed(123)
    rmidx <- sample(1:length(zdf), ceiling(length(zdf)/10), replace = FALSE)
    
    ## Convert SA1s to points, otherwise SA1s in neighbouring data supports are 
    ## also selected because they share lines, which is classified as overlapping 
    ## by over()
    SA1s_as_points <- SpatialPoints(as(SA1s, "SpatialPolygons"))
    ## Select SA1s corresponding to the selected data supports
    observed_BAU_idx <- which(!is.na(over(SA1s_as_points, zdf[rmidx, 1])))
    
    ## Do not add SA1 regions that have less than 0 families of interest (this would
    ## cause complications in model fitting)
    observed_BAU_idx <- intersect(observed_BAU_idx, 
                                  which(SA1s$Total_families_of_interest > 0))
    
    ## Combine the observed SA1s with the SA2s
    zdf <- bind(SA1s[observed_BAU_idx, ], zdf[-rmidx, ])
    
    return(zdf)
  }
}