# Source code for FRK v2 paper

This repository contains the source code for reproducing the results in "Modelling, fitting, and Prediction with Non-Gaussian Spatial and Spatio-Temporal Data using FRK". The replication may be down using a single script, `all.R`; we give detailed instructions subsequently.  


## Instructions

First, download this repo and navigate to its top-level directory within the command line (i.e., `cd` to wherever you installed the repo). 

### Data

<!---These include the Chicago crime data set and the shapefiles used in the Sydney spatial change-of-support example.--->

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

Running `all.R` populates the `results` folder with the figures and tables given in the manuscript: These can then be viewed by opening `all.html` in any web browser. To quickly establish that the code is working, low-rank versions of the models can be used: The user is prompted for their choice when running `all.R`.

<!---
"quick", low-rank versions of the models are used by default: To run the "non-quick" versions of the models and reproduce the exact results of the manuscript, set `quick = FALSE` within `all.R`. 
--->

<!---
Only a subset of the scripts need to be re-run a second time; these are clearly indicated in `all.R`, and they can be commented out to save on some computation.
--->

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

The MODIS comparison study takes a few hours to obtain the exact results of the manuscript. This was unavoidable due to the nature of the study and the necessity to give each package every opportunity to perform well. In addition to the low-rank option described above, one may easily reproduce the results from a subset of packages by editing the variable `PACKAGES` in `scripts/MODIS_control.R`. You may also wish to comment out some sections in `all.R`, as you see fit. 

<!---
### Note to Windows users

Windows users can install a Windows version of `make`: This is easy to do (see, e.g., [here](https://stackoverflow.com/a/32127632/16776594)). Once installed, the commands remain as given above. 
--->