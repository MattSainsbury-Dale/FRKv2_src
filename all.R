## Facilitates user input regardless of how this script was invoked
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

## Should we use "quick", very-low-dimensional representations of the models?
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

## Helper functions for running the scripts and providing informative output
print_start_msg <- function(section_name, section_number) {
  cat(paste("\n\n######## STARTING", section_name, "OF SECTION", section_number, "#############\n\n"))
}
  
print_run_time <- function(section_name, section_number, total_time) {
  cat(paste("\nFINISHED", section_name, "OF SECTION", section_number, "IN", round(total_time["elapsed"] / 60, 4), "MINUTES.\n"))
}
  
source_script <- function(script, section_name, section_number) {
  print_start_msg(section_name, section_number)
  total_time <- system.time(source(paste0("scripts/", script)))
  print_run_time(section_name, section_number, total_time)
  rm(list = setdiff(ls(), c("quick", "print_start_msg", "print_run_time")))
  gc()
}


## Run the scripts: 
source_script("Poisson_sim.R", "POISSON EXAMPLE", 3.1)
source_script("Negbinom_sim.R", "NEGATIVE-BINOMIAL EXAMPLE", 3.1)
source_script("Heaton.R", "HEATON COMPARISON", 3.3)
source_script("MODIS.R", "MODIS COMPARISON", 4.1)
source_script("Am.R", "AMERICIUM COMPARISON", 4.2)
source_script("Sydney.R", "SYDNEY CHANGE-OF-SUPPORT EXAMPLE", 4.3)
source_script("Chicago.R", "CHICAGO EXAMPLE", 4.4)


