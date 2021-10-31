# Source code for FRK v2 paper

This repository contains the source code for reproducing the results in "Modelling, Fitting, and Prediction with Non-Gaussian Spatial and Spatio-Temporal Data using FRK". 

## Instructions

First, download this repo and navigate to its top-level directory within the command line (i.e., `cd` to wherever you installed the repo). 

### Data

Some data sets were too large (a few hundred Mb in total) to be stored on Github: These data sets are associated with Section 4.3 (Sydney spatial change-of-support) and Section 4.4 (Chicago crime). To download them and place them into the correct folders, run `make DATA`. If you do not wish to use `make`, or if the data is not downloading as expected, please: 
- Download the  [Chicago crime data set](https://hpc.niasra.uow.edu.au/ckan/dataset/chicago_crime_dataset), name it `chicago_crime_df.Rda`, and move it to the `data` folder; and 
- download the [folders containing the Sydney SA1/SA2/SA3 region shapefiles](https://hpc.niasra.uow.edu.au/ckan/dataset/sydney_sa_regions), unzip them, and move them into `data/Sydney_shapefiles`. 

(To download the files, click "Explore" > "Go to resource".) The `data` folder should contain the following files (as well as other files automatically downloaded from Github): 

```bash
├── data
│   ├── chicago_crime_df.Rda
│   ├── Sydney_shapefiles
│   │   ├── SA1
│   │   ├── SA2
│   │   ├── SA3
```

Note that we have included checks at the beginning of the replication script to ensure that the user is immediately notified if these files are not present.

### Dependencies

When running the replication script, the user will be asked if they wish to install package dependencies. If they choose to do so, they will then be asked if pre-existing packages should be re-installed with the correct version numbers as given in `dependencies.txt` (this option is only recommended for use if there is a problem with the latest version of the packages).

### Results

The replication script is `run_all.sh`, invoked using `bash run_all.sh`. Alternatively, Windows users may use `run_all.bat`. The replication script populates the `results` folder with the figures and tables given in the manuscript: These can then be viewed by opening `results_all.html` in any web browser. To quickly establish that the code is working, very-low-dimensional representations of the models can be used: The user is prompted for their choice when running the replication script. Note that some of the results will look very different when using these low-dimensional representations (particularly the results for Section 4.1). 

Please note that the replication material for Section 3.3, Table 4, reproduces the results for row FRK v2 only; the results for other rows are taken from Table 3 in Heaton et al. (2019). 

#### Hardware requirements

You will need at least 32GB of RAM (or RAM + swap) to run the very-low-dimensional representations of the models, and you will need at least 64GB of RAM (or RAM + swap) to run the full models. 

#### Run times

The following are estimates of the expected run time needed to reproduce the full results of the manuscript: 

- Section 3.1: ~5 minutes
- Section 3.2: ~5 minutes
- Section 3.3: ~30 minutes
- Section 4.1: ~15 hours
- Section 4.2: ~5 minutes
- Section 4.3: ~5 minutes
- Section 4.4: ~30 minutes

Please note that the long run-time for Section 4.1 (MODIS comparison study) was unavoidable due to the nature of the study and the necessity to provide each package with every opportunity to perform well. In addition to the option of using low-dimensional representations as described above, one may easily reproduce the results from a subset of packages by editing the variable `PACKAGES` in `scripts/MODIS.R`. We thank reviewers for their patience and understanding. 

Note that the replication script is clearly presented and commented; hence, one may easily "comment out" sections to produce a subset of the results. (Comments in `.sh` files are made with `##`, and comments in `.bat` files are made using `::`.)


#### Real-world example

For users wishing to "try out" the package on a real-world example, we suggest the Sydney spatial change-of-support example presented in Section 4.3.  

<!---
### Note to Windows users

Windows users can install a Windows version of `make`: This is easy to do (see, e.g., [here](https://stackoverflow.com/a/32127632/16776594)). Once installed, the commands remain as given above. 
--->