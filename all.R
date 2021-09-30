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

## Should we use "quick", very-low-dimensional representations of the models?
quick <- user_decision("Do you wish to use very-low-dimensional representations of the models to quickly establish that the code is working? Note that the generated results and plots will not exactly match those in the manuscript if you reply 'y', and you will need at least 32GB of RAM and/or swap space if you reply 'n'. (y/n)")
quick <- quick == "y" # Convert to Boolean

if (!quick) {
  Chicago_nres <- user_decision("Do you want to use nres = 2 or nres = 3 in the Chicago section (Section 4.4 in the paper)? The paper used nres = 3, but it is computationally demanding and requires at least 64 GB of RAM, and similar results can be obtained with nres = 2. (Please enter 2 or 3)", 
                                allowed_answers = c(2, 3))
  Chicago_nres <- as.numeric(Chicago_nres)
  
  } else {
  Chicago_nres <- 1
}

if (Chicago_nres == 3) {
  cat("You have chosen to use nres = 3 in the Chicago section. Depending on your hardware, the Chicago section may need a new R session to avoid memory problems: If you encounter memory problems during this section, please comment out all sections bar the Chicago section, and then re-run all.R.\n")
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

## Helper functions for running the scripts and providing informative output
print_start_msg <- function(section_name, section_number) {
  cat(paste("\n\n######## STARTING", section_name, "OF SECTION", section_number, "#############\n\n"))
}
  
print_run_time <- function(section_name, section_number, total_time) {
  cat(paste("\nFinished", tolower(section_name), "of Section", section_number, "in", round(total_time["elapsed"] / 60, 4), "minutes.\n"))
}
  
source_wrapper <- function(script, section_name, section_number) {
  print_start_msg(section_name, section_number)
  total_time <- system.time(source(paste0("scripts/", script)))
  print_run_time(section_name, section_number, total_time)
  rm(list = setdiff(ls(), c("quick", "print_start_msg", "print_run_time")))
  invisible(gc())
}


## Run the scripts: 
source_wrapper("Poisson_sim.R", "POISSON EXAMPLE", 3.1)
source_wrapper("Negbinom_sim.R", "NEGATIVE-BINOMIAL EXAMPLE", 3.2)
source_wrapper("Heaton.R", "HEATON COMPARISON", 3.3)
source_wrapper("MODIS.R", "MODIS COMPARISON", 4.1)
source_wrapper("Am.R", "AMERICIUM COMPARISON", 4.2)
source_wrapper("Sydney.R", "SYDNEY CHANGE-OF-SUPPORT EXAMPLE", 4.3)
source_wrapper("Chicago.R", "CHICAGO EXAMPLE", 4.4)


