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
    tmp <- paste(allowed_answers, collapse = " or ")
    cat(paste0("Please enter ", tmp, ".\n"))
    answer <- user_decision(prompt, allowed_answers = allowed_answers)
  }
  
  return(answer)
}

## Install dependencies:
install_depends <- user_decision("Do you want to automatically install package dependencies? (y/n)")
if (install_depends == "y") {
  install_exact_versions <- user_decision("Do you want to ensure that all package versions are as given in dependencies.txt? (y/n)")
  install_exact_versions <- install_exact_versions == "y" # Convert to Boolean
  
  if (install_exact_versions) {
    cat("When changing the packages to the versions specified in dependencies.txt, please use your discretion when answering the question “Which would you like to update?”.  Updating all packages (i.e., option 3) may cause errors.")
  }
  
  source("scripts/Dependencies_install.R")
}


## Check that the data has been downloaded: 
downloaded_correctly <- c(
  file.exists("./data/chicago_crime_df.Rda"),
  # Check that the Sydney SA directories exist:
  file.exists("./data/Sydney_shapefiles/SA1"), 
  file.exists("./data/Sydney_shapefiles/SA2"), 
  file.exists("./data/Sydney_shapefiles/SA3"),
  # Check that the actual shapefiles exist:
  file.exists("./data/Sydney_shapefiles/SA1/SA1_2011_AUST.shp"), 
  file.exists("./data/Sydney_shapefiles/SA2/SA2_2011_AUST.shp"), 
  file.exists("./data/Sydney_shapefiles/SA3/SA3_2011_AUST.shp")
)
names(downloaded_correctly) <- c("SA1", "SA2", "SA3", "chicago_crime_df.Rda")
if(!all(downloaded_correctly)) {
  stop(paste0(
    "Some files have not been downloaded or are in the wrong place.\nProblematic files: \n", 
    names(downloaded_correctly)[!downloaded_correctly], 
    "\nIf on Linux/Unix, run 'make DATA', otherwise please see the README for download instructions."
  ))
}
