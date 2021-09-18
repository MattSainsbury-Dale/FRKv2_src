CRANMIRROR <- "https://cran.csiro.au" 

## Find the packages used throughout this repo using the package renv
if (!("renv" %in% rownames(installed.packages()))) 
  install.packages("renv", repos = CRANMIRROR)
depends <- renv::dependencies()
depends <- unique(depends$Package)

## Assuming all of the packages are installed, find the package versions that 
## present on THIS system: This was mainly done so that we have a record of the 
## package versions at the time of the manuscript submission, and hopefully 
## we can ensure reproducible results in the future. 
pkg_versions <- sapply(depends, function(pkg) as.character(packageVersion(pkg)))

## Save the list of dependencies to a text file
write.table(data.frame(depends, pkg_versions), "dependencies.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)




