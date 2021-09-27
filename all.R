# Facilitates user input regardless of how this script was invoked
user_decision <- function(prompt) {
  
  if (interactive()) {
    answer <- readline(prompt)
  } else {
    cat(prompt)
    answer <- readLines("stdin", n = 1)
  }
  
  answer <- tolower(answer)
  
  if (answer != "y" && answer != "n") {
    cat("Please enter y or n.\n")
    answer <- user_decision(prompt)
  }
  
  return(answer)
}

## Should we use "quick", low-rank versions of the models?
quick <- user_decision("Do you wish to use very-low-dimensional representations of the models to quickly establish that the code is working (note that the generated results and plots will not exactly match those in the manuscript if you reply 'y', and you will need at least 32GB of RAM and/or swap space if you reply 'n')? (y/n)\n")
quick <- quick == "y" # Convert to Boolean


## Install dependencies:
install_depends <- user_decision("Do you want to automatically install package dependencies? (y/n)\n")
if (install_depends == "y") {
  install_exact_versions <- user_decision("Do you want to ensure that all package versions are as given in dependencies.txt? (y/n)\n")
  install_exact_versions <- install_exact_versions == "y" # Convert to Boolean
  source("scripts/Dependencies_install.R")
}


## Check that the data has been downloaded: 
downloaded_correctly <- c(
  file.exists("./data/Sydney_shapefiles/SA1/"), 
  file.exists("./data/Sydney_shapefiles/SA2/"), 
  file.exists("./data/Sydney_shapefiles/SA3/"), 
  file.exists("./data/chicago_crime_df.Rda")
  )
names(downloaded_correctly) <- c("SA1", "SA2", "SA3", "chicago_crime_df.Rda")
if(!all(downloaded_correctly)) {
  stop(paste0(
    "Some files have not been downloaded or are in the wrong place.\nProblematic files: \n", 
    names(downloaded_correctly)[!downloaded_correctly], 
    "\nIf on Linux/Unix, run 'make DATA', otherwise please see the README for download instructions."
    ))
}




cat("\n\n######## STARTING POISSON EXAMPLE OF SECTION 3.1 #############\n\n")
source("scripts/Poisson_sim.R")
rm(list = setdiff(ls(), "quick"))
gc()

cat("\n\n######## STARTING NEGATIVE-BINOMIAL EXAMPLE OF SECTION 3.2 #############\n\n")
source("scripts/Negbinom_sim.R")
rm(list = setdiff(ls(), "quick"))
gc()

cat("\n\n######## STARTING HEATON COMPARISON OF SECTION 3.3 #############\n\n")
source("scripts/Heaton.R")
rm(list = setdiff(ls(), "quick"))
gc()

cat("\n\n######## STARTING MODIS COMPARISON OF SECTION 4.1 #############\n\n")
source("scripts/MODIS.R")
rm(list = setdiff(ls(), "quick"))
gc()

cat("\n\n######## STARTING AMERICIUM COMPARISON OF SECTION 4.2 #############\n\n")
source("scripts/Am.R")
rm(list = setdiff(ls(), "quick"))
gc()

cat("\n\n######## STARTING SYDNEY EXAMPLE OF SECTION 4.3 #############\n\n")
source("scripts/Sydney.R")
rm(list = setdiff(ls(), "quick"))
gc()

cat("\n\n######## STARTING CHICAGO EXAMPLE OF SECTION 4.4 #############\n\n")
source("scripts/Chicago.R")
rm(list = setdiff(ls(), "quick"))
gc()
