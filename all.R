# Code that is unaffected by the variable "quick":
source("scripts/Poisson_sim.R")
source("scripts/Sydney.R")
source("scripts/Negbinom_sim.R")
source("scripts/Am_GZ.R")
source("scripts/Am_data.R")
source("scripts/Am_blocks.R")
source("scripts/Am_BAUs.R")
source("scripts/Am_data_blocks_plot.R")
source("scripts/Am_georob.R")
source("scripts/Am_FRK.R")
source("scripts/Am_comparison_plot.R")

# Set the variable "quick":
quick <- TRUE

# Code that is affected by the variable "quick":
source("scripts/Heaton.R")
source("scripts/MODIS_control.R")
source("scripts/Chicago.R") 