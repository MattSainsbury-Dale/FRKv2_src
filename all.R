## Set the variable "quick":
quick <- TRUE

## Install dependencies:
## (See here (https://stackoverflow.com/a/44013192/16776594) for a trick to make 
## tprompting work from the command line)

if (interactive()) {
  answer1 <- readline(cat("Do you want to automatically install package dependencies? (Y/N)\n"))
  answer1 <- toupper(answer1)
  
  if (answer1 != "Y" && answer1 != "N") {
    stop("Please enter Y or N.")
  } else if (answer1 == "Y") {
    answer2 <- readline(cat("Do you want to ensure that all package versions are as given in dependencies.txt? (Y/N)\n"))
    answer2 <- toupper(answer2)
    if (answer2 != "Y" && answer2 != "N") {
      stop("Please enter Y or N.")
    } else if (answer == "N") {
      install_correct_versions <- FALSE
    } else if (answer == "Y") {
      install_correct_versions <- TRUE
    }
    source("scripts/Dependencies_install.R")
  }
} else {
  ## Just install the packages but do not enforce the correct versions 
  install_correct_versions <- FALSE
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
    "\nPlease see the README for download instructions."
    ))
}


cat("\n\n######## STARTING POISSON EXAMPLE OF SECTION 3.1 #############\n\n")
source("scripts/Poisson_sim.R")

cat("\n\n######## STARTING NEGATIVE-BINOMIAL EXAMPLE OF SECTION 3.2 #############\n\n")
source("scripts/Negbinom_sim.R")

cat("\n\n######## STARTING HEATON COMPARISON OF SECTION 3.3 #############\n\n")
source("scripts/Heaton.R")

cat("\n\n######## STARTING MODIS COMPARISON OF SECTION 4.1 #############\n\n")
source("scripts/MODIS_control.R")

cat("\n\n######## STARTING AMERICIUM COMPARISON OF SECTION 4.2 #############\n\n")
source("scripts/Am_GZ.R")
source("scripts/Am_data.R")
source("scripts/Am_blocks.R")
source("scripts/Am_BAUs.R")
source("scripts/Am_data_blocks_plot.R")
source("scripts/Am_georob.R")
source("scripts/Am_FRK.R")
source("scripts/Am_comparison_plot.R")

cat("\n\n######## STARTING SYDNEY EXAMPLE OF SECTION 4.3 #############\n\n")
source("scripts/Sydney.R")

cat("\n\n######## STARTING CHICAGO EXAMPLE OF SECTION 4.4 #############\n\n")
source("scripts/Chicago.R") 


