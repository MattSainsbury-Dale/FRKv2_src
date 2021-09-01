# ---- Install dependencies ----


warning(
"This script does not ensure that the correct package version is installed. 
This should not be a problem, as R is a stable language and most packages are 
backwards compatible. However, if problems with the code occur, please manually 
install the package versions as used at the time of writing (see dependencies.txt)."
)

depends_and_pkg_versions <- read.table("dependencies.txt", header = FALSE)
depends <- as.character(depends_and_pkg_versions[, 1])

## Find the packages not yet installed
new_packages <- depends[!(depends %in% rownames(installed.packages()))] 

## If FRK is already present but not FRK v.2, 
## add it to the list of packages to be installed. 
if (!("FRK" %in% new_packages) && packageVersion("FRK") < 2) 
  new_packages <- c(new_packages, "FRK")

if(length(new_packages)) {
  warning("Uninstalled package dependencies are being installed automatically 
           via the script src/Dependencies_install.R")
  install.packages(new_packages, repos = CRANMIRROR)
} else {
  message("\nNo package dependecies to be installed.\n")
}