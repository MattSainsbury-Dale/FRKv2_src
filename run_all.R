## Data and dependencies
source("setup.R")

## Should we use "quick", very-low-dimensional representations of the models?
quick <- user_decision("Do you wish to use very-low-dimensional representations of the models to quickly establish that the code is working? Note that the generated results and plots will not exactly match those in the manuscript if you reply 'y', and you will need at least 64GB of RAM and/or swap space if you reply 'n'. (y/n)")
quick <- quick == "y" # Convert to Boolean

if (!quick) {
  Chicago_nres <- user_decision("Do you want to use nres = 2 or nres = 3 in the Chicago section (Section 4.4 in the paper)? The paper used nres = 3, but it is computationally demanding, and similar results can be obtained with nres = 2. (Please enter 2 or 3)", 
                                allowed_answers = c(2, 3))
  Chicago_nres <- as.numeric(Chicago_nres)
  
} else {
  Chicago_nres <- 1
}

if (Chicago_nres == 3) {
  cat("You have chosen nres = 3... Depending on your hardware, the Chicago section may need a new R session to avoid memory problems: If you encounter memory problems during this section, please comment out all sections bar the Chicago section, and then re-run all.R.\n")
}

source("scripts/Utility_fns.R")


## Run the scripts: 
source_wrapper("Poisson_sim.R", "POISSON EXAMPLE", 3.1)
source_wrapper("Negbinom_sim.R", "NEGATIVE-BINOMIAL EXAMPLE", 3.2)
source_wrapper("Heaton.R", "HEATON COMPARISON", 3.3)
source_wrapper("MODIS.R", "MODIS COMPARISON", 4.1)
source_wrapper("Am.R", "AMERICIUM COMPARISON", 4.2)
source_wrapper("Sydney.R", "SYDNEY CHANGE-OF-SUPPORT EXAMPLE", 4.3)
source_wrapper("Chicago.R", "CHICAGO EXAMPLE", 4.4)


