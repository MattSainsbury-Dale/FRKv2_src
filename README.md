# Source code for FRK v.2 manuscript

This repository contains the source code to replicate the figures and results of the FRK v.2 manuscript. 

## Instructions

To reproduce the results please download this repository (see [here](https://superuser.com/a/1309684) for steps to download a repository). In terminal, navigate to the top directory of the repository. It is important that the working directory is the top directory because the scripts use references to ./src/, ./img/, ./results/, and so on. Then, one may use the makefile to automatically populate the img/ and results/ directories. The targets in the makefile correspond to subsections of Sections 3 and 4 of the manuscript. Specially, the main targets are
- `all`	Produces all of the figures and results in the manuscript
  - `Poisson_sim` Produces the figures and results of Section 3.1 
  - `Negbinom_sim` Produces the figures and results of Section 3.2
  - `Heaton` Produces the FRK v.2 entry of the Heaton et al. (2019) comparison study table shown in Section 3.4
  - `MODIS` Produces the figures and results of Section 4.1 (see below for some considerations)
  - `Am` Produces the figures and results of Section 4.2
  - `Sydney` Produces the figures and results of Section 4.3
  - `Chicago` Produces the figures and results of Section 4.4
To wipe the populated directories after a call to make, simply enter `make clean`.

### Data

Some of the data files were too large to be stored on Github. Before calling `make`, please run the following commands in terminal (assuming you have navigated to the top directory of this repo) to download these files and place them in the appropriate place.
- `curl https://dl.dropboxusercontent.com/s/cwlziqavrpe7ode/chicago_crime_df -O data/chicago_crime_df.Rda`
- `curl https://dl.dropboxusercontent.com/s/v75fpk0r6pig8o4/SA1_2011_AUST.shp -O data/Sydney_shapefiles/SA1/SA1_2011_AUST.shp`
- `curl https://dl.dropboxusercontent.com/s/m5ocngvauxosz50/SA2_2011_AUST.shp -O data/Sydney_shapefiles/SA2/SA2_2011_AUST.shp`
- `curl https://dl.dropboxusercontent.com/s/m4qigm9s4ijn4bj/SA3_2011_AUST.shp -O data/Sydney_shapefiles/SA3/SA3_2011_AUST.shp`




### MODIS

The MODIS study is computationally demanding (upwards of 10 hours). To simply test that the code works, one may use extremely low rank versions of the models by setting LOWRANK = TRUE within src/MODIS_control.R. In fact, LOWRANK = TRUE is the default, so that people do not inadvertently get stuck in some lengthy computations. Further, one can reproduce the results from a subset of packages by editing the object PACKAGES in src/MODIS_control.R.

#### References

* Heaton MJ, Datta A, Finley AO, Furrer R, Guinness J, Guhaniyogi R, Gerber F, Gramacy
RB, Hammerling D, Katzfuss M, Lindgren F, Nychka DW, Sun F, Zammit-Mangion A
(2019). “A case study competition among methods for analyzing large spatial data.” Journal
of Agricultural, Biological and Environmental Statistics, 24, 398–425.