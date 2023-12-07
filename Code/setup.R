# ---- Helper functions for this script ----

## Facilitates user input regardless of how this script was invoked
user_decision <- function(prompt, allowed_answers = c("y", "n")) {
  
  if (interactive()) {
    answer <- readline(paste0(prompt, "\n"))
  } else {
    cat(prompt)
    answer <- readLines("stdin", n = 1)
  }
  
  answer <- tolower(answer)
  
  if (!(answer %in% allowed_answers)) {
    cat(paste0("Please enter ", paste(allowed_answers, collapse = " or "), ".\n"))
    answer <- user_decision(prompt, allowed_answers = allowed_answers)
  }
  
  return(answer)
}

# ---- Install dependencies ----

install_depends <- user_decision("Do you want to automatically install package dependencies? (y/n)")
if (install_depends == "y") {
  install.packages("dggrids", repos="https://andrewzm.github.io/dggrids-repo", type = "source")
  install.packages("INLA", repos="https://inla.r-inla-download.org/R/stable")
  install.packages(c("maps", "devtools", "dplyr", "DHARMa", "foreach", "FRK", "fmesher", "ggmap", "ggplot2", "ggpubr", "gstat", "htmltab", "htmlTable", "mgcv" , "plyr", "pROC" , "R.utils" , "raster", "renv" , "reshape2", "scales", "scoringRules", "sp", "sf", "spacetime", "spBayes", "splancs" , "spNNGP", "stringr", "tibble", "tidyr"), 
                   repos = "https://cloud.r-project.org")
}


# ---- Check that the data has been downloaded correctly ----

downloaded_correctly <- c(
  file.exists("data/chicago_crime_df.Rda"),
  # Check that the Sydney SA directories exist:
  file.exists("data/Sydney_shapefiles/SA1"), 
  file.exists("data/Sydney_shapefiles/SA2"), 
  file.exists("data/Sydney_shapefiles/SA3"),
  # Check that the shapefiles exist:
  file.exists("data/Sydney_shapefiles/SA1/SA1_2011_AUST.shp"), 
  file.exists("data/Sydney_shapefiles/SA2/SA2_2011_AUST.shp"), 
  file.exists("data/Sydney_shapefiles/SA3/SA3_2011_AUST.shp")
)
names(downloaded_correctly) <- c("SA1", "SA2", "SA3", "chicago_crime_df.Rda")
if(!all(downloaded_correctly)) {
  stop(paste0(
    "Some files have not been downloaded or are in the wrong place.\nProblematic files: \n", 
    names(downloaded_correctly)[!downloaded_correctly], 
    "\nIf on a unix system, run 'bash DATA.sh', otherwise please see the README for download instructions."
  ))
}
