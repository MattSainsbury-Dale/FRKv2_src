# Source code for FRK v2 paper

This repository contains the source code to replicate the figures and results of the FRK v2 paper. 

## Instructions


<!---
First download this repo and navigate to its *top-level directory* within terminal (i.e., `cd` to wherever you installed the repo). Then, one may use `make` to automatically run the source code in the `scripts` folder, populating the `img` and `results` directories with the figures and results of the manuscript. The main targets in the `makefile` correspond to subsections of Sections 3 and 4 of the manuscript:
- `make all`	Produces all of the figures and results in the manuscript
  - `make Poisson_sim` Produces the figures and results of Section 3.1 
  - `make Negbinom_sim` Produces the figures and results of Section 3.2
  - `make Heaton` Produces the FRK v2 entry for the comparison study of Heaton et al. (2019), shown in Section 3.4
  - `make MODIS` Produces the figures and results of Section 4.1 (see below for considerations regarding run time)
  - `make Am` Produces the figures and results of Section 4.2
  - `make Sydney` Produces the figures and results of Section 4.3 (see below for considerations regarding data)
  - `make Chicago` Produces the figures and results of Section 4.4 (see below for considerations regarding data and run time)
--->

  
To wipe the populated directories, run `make clean` (or manually delete the contents of `img/` and `results/`).

### Dependencies

The file `dependencies.txt` contains the `R` package dependencies with version numbers at the time that the manuscript was prepared. To automatically install the required packages that are not already present on the user's system, one may run `make INSTALL_DEPENDS` or run `scripts/Dependencies_install.R`. Note that this uses the CRAN mirror run by the CSIRO in Australia.

### Data

Some data sets used in this analysis are too large to be stored on Github: These include the Chicago crime data set and the shapefiles used in the Sydney spatial change-of-support example. These data sets are instead conveniently stored [here](https://hpc.niasra.uow.edu.au/ckan/ar/dataset/chicago_crime_dataset) and [here](https://hpc.niasra.uow.edu.au/ckan/ar/dataset/sydney_sa_regions): To download them and place them into the correct folders, run `make DATA`. FIXME: WHAT ABOUT WINDOWS USERS?

### A note on long run times

The MODIS study is computationally demanding, taking upwards of 10 hours to reproduce the results of the manuscript. To simply test that the code works, one may use extremely low-rank versions of the models by setting `LOWRANK = TRUE` within `scripts/MODIS_control.R`. Further, one can reproduce the results from a subset of packages by editing the object `PACKAGES` in `scripts/MODIS_control.R`. The Chicago application study is also computationally demanding: A low-rank version of the model can be used by setting `nres = 1` within `scripts/Chicago.R`.


### Note to Windows users

Windows user may need to install a Windows version of `make`: This is easy to do (see, e.g., [here](https://stackoverflow.com/a/32127632/16776594)). Once installed, the commands remain as given above. 