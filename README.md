# Source code for FRK v2 paper

This repository contains the source code for reproducing the results in "Modelling, Fitting, and Prediction with Non-Gaussian Spatial and Spatio-Temporal Data using FRK". 

## Instructions

First, download this repo and navigate to its top-level directory within the command line (i.e., `cd` to wherever you installed the repo). 

### Data

Some data sets were too large to be stored on Github (a few hundred Mb in total): To download them and place them into the correct folders, run `make DATA`. If you do not wish to use `make`, or if the data is not downloading as expected, please: 
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

When running the replication script, the user will be asked if they wish to install package dependencies; if they choose to do so, they will then be asked if pre-existing packages should be re-installed with the correct version numbers as given in `dependencies.txt`. 

### Results

The replication script is `run_all.sh`, invoked using `bash run_all.sh`. Alternatively, Windows users may use `run_all.bat`. The replication script populates the `results` folder with the figures and tables given in the manuscript: These can then be viewed by opening `results_all.html` in any web browser. To quickly establish that the code is working, very-low-dimensional representations of the models can be used: The user is prompted for their choice when running the replication script. Note that some of the results will look very different when using these low-dimensional representations (particularly the results for Section 4.1).

#### Hardware requirements

You will need at least 32GB of RAM (or RAM + swap) to run the very-low-dimensional representations of the models, and you will need at least 64GB of RAM (or RAM + swap) to run the full models. 

#### A note on long run times

The MODIS comparison study takes several hours to obtain the exact results of the manuscript. This was unavoidable due to the nature of the study and the necessity to give each package every opportunity to perform well. In addition to the option of using low-dimensional representations as described above, one may easily reproduce the results from a subset of packages by editing the variable `PACKAGES` in `scripts/MODIS.R`. You may also wish to comment out some sections in the replication script, as you see fit. 

<!---
### Note to Windows users

Windows users can install a Windows version of `make`: This is easy to do (see, e.g., [here](https://stackoverflow.com/a/32127632/16776594)). Once installed, the commands remain as given above. 
--->