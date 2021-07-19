# Source code for FRK v.2 manuscript

This repository contains the source code to replicate the figures and results of the FRK v.2 manuscript. 

## Instructions

To reproduce the results please download this repository (see [here](https://superuser.com/a/1309684) for steps to download a repository). In terminal, navigate to the top directory of the repository. It is important that the working directory is the top directory because the scripts use references to src/, img/, results/, and so on. Then, one may use `make` to automatically run the source code and populate the img/ and results/ directories with the figures and results of the manuscript. The targets in the makefile correspond to subsections of Sections 3 and 4 of the manuscript. 
- `make all`	Produces all of the figures and results in the manuscript
  - `make Poisson_sim` Produces the figures and results of Section 3.1 
  - `make Negbinom_sim` Produces the figures and results of Section 3.2
  - `make Heaton` Produces the FRK v.2 entry of the Heaton et al. (2019) comparison study table shown in Section 3.4
  - `make MODIS` Produces the figures and results of Section 4.1 (see below for some considerations)
  - `make Am` Produces the figures and results of Section 4.2
  - `make Sydney` Produces the figures and results of Section 4.3
  - `make Chicago` Produces the figures and results of Section 4.4
  
To wipe the populated directories, enter `make clean`.

### Dependencies

The file DEPENDS.txt contains the `R` package dependencies for this repo. To automatically install the required packages (if not already present), please run `make DEPENDS`, which calls the `R` script src/DEPENDS.R. Note that this script uses the CRAN mirror run by the CSIRO (in Australia), so you may want to change this.

### Data

Some data sets used in this analysis are too large to be stored on Github: These include the Chicago crime data set and the shapefiles used in the Sydney spatial change-of-support example. These are instead stored on the NIASRA HPC: To download and place them into the data/ sub-directory, please run `make DATA` within terminal (from the top level of the directory). 

### Long run times

The MODIS study is computationally demanding, taking upwards of 10 hours to reproduce the results of the manuscript. To simply test that the code works, one may use extremely low-rank versions of the models by setting `LOWRANK = TRUE` within src/MODIS_control.R. Further, one can reproduce the results from a subset of packages by editing the object `PACKAGES` in src/MODIS_control.R. The Chicago application study is also computationally demanding: A low-rank version of the model can be used by setting nres = 1 instead of nres = 3 within src/Chicago.R.
