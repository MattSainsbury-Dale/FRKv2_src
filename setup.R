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

