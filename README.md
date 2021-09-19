# Source code for FRK v2 paper

This repository contains the source code to replicate the results of the FRK v2 paper. 

## Instructions

First, download this repo and navigate to its top-level directory within the command line (i.e., `cd` to wherever you installed the repo). 

#### Dependencies


The file `dependencies.txt` contains the package dependencies and version numbers. To automatically install the dependencies that are not already present on your system, enter `Rscript scripts/Dependencies_install.R`. Note that this uses the CRAN mirror hosted by the CSIRO in Australia. 

<!---run `make INSTALL_DEPENDS` or --->

#### Data

<!---These include the Chicago crime data set and the shapefiles used in the Sydney spatial change-of-support example.--->

Some data sets were too large to be stored on Github: To download them and place them into the correct folders, run `make DATA`. If you do not wish to use `make`, please: 
- Download the  [Chicago crime data set](https://hpc.niasra.uow.edu.au/ckan/ar/dataset/chicago_crime_dataset), name it `chicago_crime_df.Rda`, and move it to the `data` folder; and 
- download the [Sydney shapefiles](https://hpc.niasra.uow.edu.au/ckan/ar/dataset/sydney_sa_regions), unzip them, and move them into `data/Sydney_shapefiles`.

#### Results

To reproduce the results of the manuscript, simply run `Rscript all.R`. This populates the `results` folder with the figures and tables given in the manuscript; these can then be viewed by opening `all.html` in any web browser. To alleviate the computational burden of reviewing this code, "quick", low-rank versions of the models are used by default: This is done to quickly establish that the code is in working order. To run the "non-quick" versions of the models and reporduce the exact results of the manuscript, set `quick = FALSE` within `all.R`. Only a subset of the scripts need to be re-run a second time; these are clearly indicated in `all.R`, and they can be commented out to save on some computation.

<!---
To alleviate long run-time issues, we have provided an option to use low-rank versions of the models: This is done by manually setting `quick = TRUE` within `all.R`. Our envisioned workflow for a reviewer is to first run the code with `quick = TRUE` to quickly establish that it is in working order, and then re-run it with `quick = FALSE`. Only a subset of the scripts need to be re-run a second time; these are clearly indicated in `all.R`, and they can be commented out to save on some computation.
--->

<!---
We provide two convenient options for reproducing the results of the manuscript: 

- Run `all.R` 
- Run `make all` 

Both options populate the `results` folder with the figures and tables used in the manuscript; these can then be viewed by opening `all.html` in any web browser. To alleviate long run-time issues, we have provided an option to use low-rank versions of the models: This is done by manually setting `quick=TRUE` within `all.R`, or by running `make all quick=TRUE`. Our envisioned workflow for a reviewer is to first run the code with `quick=TRUE` to quickly establish that it is in working order, and then re-run it with `quick=FALSE`. If using `make`, only a subset of the scripts will be re-run a second time; if you are manually running `all.R`, you can comment-out the scripts that do not depend on the `quick` variable (these are clearly indicated in the script).
--->


<!---
First download this repo and navigate to its *top-level directory* within terminal (i.e., `cd` to wherever you installed the repo). Then, one may use `make` to automatically run the source code in the `scripts` folder, populating the `img` and `results` directories with the figures and results of the manuscript. The main benefit of using `make` is that it tracks the timestamps of created files, which can help to avoid unnecessary computation. The main targets in the `makefile` correspond to subsections of Sections 3 and 4 of the manuscript:
- `make all`	Produces all of the figures and results in the manuscript
  - `make Poisson_sim` Produces the figures and results of Section 3.1 
  - `make Negbinom_sim` Produces the figures and results of Section 3.2
  - `make Heaton` Produces the FRK v2 entry for the comparison study of Heaton et al. (2019), shown in Section 3.4
  - `make MODIS` Produces the figures and results of Section 4.1 (see below for considerations regarding run time)
  - `make Am` Produces the figures and results of Section 4.2
  - `make Sydney` Produces the figures and results of Section 4.3 (see below for considerations regarding data)
  - `make Chicago` Produces the figures and results of Section 4.4 (see below for considerations regarding data and run time)
--->

  
<!---To wipe the populated directories, run `make clean` or manually delete the contents of `results/`.--->



#### A note on long run times

The MODIS study takes a few hours if `quick = FALSE`. This was unavoidable due to the nature of the study and the necessity to give each package every opportunity to perform well. In addition to the `quick = TRUE` option described above, one may easily reproduce the results from only a subset of packages by editing the variable `PACKAGES` in `scripts/MODIS_control.R`. 

<!---
### Note to Windows users

Windows users can install a Windows version of `make`: This is easy to do (see, e.g., [here](https://stackoverflow.com/a/32127632/16776594)). Once installed, the commands remain as given above. 
--->