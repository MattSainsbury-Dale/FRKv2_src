# Source code for FRK v2 paper

This repository contains the source code for reproducing the results in "Modelling, fitting, and Prediction with Non-Gaussian Spatial and Spatio-Temporal Data using FRK". The replication may be down using a single script, `all.R`; we give detailed instructions subsequently.  

## Instructions

First, download this repo and navigate to its top-level directory within the command line (i.e., `cd` to wherever you installed the repo). 

### Data

Some data sets were too large to be stored on Github (a few hundred Mb in total): To download them and place them into the correct folders, run `make DATA`. If you do not wish to use `make`, or if the data is not downloading as expected, please: 
- Download the  [Chicago crime data set](https://hpc.niasra.uow.edu.au/ckan/dataset/chicago_crime_dataset), name it `chicago_crime_df.Rda`, and move it to the `data` folder; and 
- download the [folders containing the Sydney SA1/SA2/SA3 region shapefiles](https://hpc.niasra.uow.edu.au/ckan/dataset/sydney_sa_regions), unzip them, and move them into `data/Sydney_shapefiles`. 

(To download the files, click "Explore" > "Go to resource")

The `./data/` folder should contain the following files (as well as other files that come with the repo), and we have included checks in `all.R` to ensure that this is the case:

```bash
├── data
│   ├── chicago_crime_df.Rda
│   ├── Sydney_shapefiles
│   │   ├── SA1
│   │   ├── SA2
│   │   ├── SA3
```

### Dependencies

When running `all.R`, the user will be asked if they wish to install package dependencies; if they choose to do so, they will then be asked if pre-existing packages should be re-installed with the correct version numbers as given in `dependencies.txt`. 

### Results

Running `all.R` (e.g., with `Rscript all.R`) populates the `results` folder with the figures and tables given in the manuscript: These can then be viewed by opening `all.html` in any web browser. To quickly establish that the code is working, very-low-dimensional representations of the models can be used: The user is prompted for their choice when running `all.R`. Note that you will need at least 32GB of RAM and/or swap space to run the full models. 

#### A note on long run times

The MODIS comparison study takes a few hours to obtain the exact results of the manuscript. This was unavoidable due to the nature of the study and the necessity to give each package every opportunity to perform well. In addition to the low-rank option described above, one may easily reproduce the results from a subset of packages by editing the variable `PACKAGES` in `scripts/MODIS.R`. You may also wish to comment out some sections in `all.R`, as you see fit. 

<!---
### Note to Windows users

Windows users can install a Windows version of `make`: This is easy to do (see, e.g., [here](https://stackoverflow.com/a/32127632/16776594)). Once installed, the commands remain as given above. 
--->