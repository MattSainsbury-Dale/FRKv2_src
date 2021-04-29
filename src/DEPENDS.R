CRANMIRROR <- "https://cran.csiro.au" 

## Find the packages used throughout this repo using the package renv
if (!("renv" %in% rownames(installed.packages()))) 
  install.packages("renv", repos = CRANMIRROR)
DEPENDS <- renv::dependencies()
DEPENDS <- unique(DEPENDS$Package)

## Save the list of dependencies to a text file
write.table(DEPENDS, "DEPENDS.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

## Find the packages not yet installed
new_packages <- DEPENDS[!(DEPENDS %in% rownames(installed.packages()))] 

## If FRK is already present but not FRK v.2, 
## add it to the list of packages to be installed. 
if (!("FRK" %in% new_packages) && packageVersion("FRK") < 2)
  new_packages <- c(new_packages, "FRK")

if(length(new_packages)) {
  warning("Uninstalled package dependencies are being installed automatically by the script src/DEPENDS.R")
  install.packages(new_packages, repos = CRANMIRROR)
}