# ---- Install dependencies ----

if (!("devtools" %in% rownames(installed.packages()))) {
  install.packages("devtools", repos = CRANMIRROR)
}

CRANMIRROR <- "https://cran.csiro.au" 

tmp <- read.table("dependencies.txt", header = FALSE)
pkg_versions <- setNames(as.character(tmp[, 2]), tmp[, 1])
rm(tmp)

## Find the packages not yet installed and add them to the list
installed_idx <- names(pkg_versions) %in% rownames(installed.packages())
new_packages  <- names(pkg_versions)[!(installed_idx)] 

if (exists("install_correct_versions") && install_correct_versions) {
  ## Find the packages that are installed, but not the correct version
  installed_pkg_versions <- sapply(names(pkg_versions)[installed_idx], function(pkg) as.character(packageVersion(pkg)))
  idx          <- installed_pkg_versions != pkg_versions[installed_idx]
  already_installed_pkgs_different_versions <- names(installed_pkg_versions)[idx]
  # new_packages <- c(new_packages, names(installed_pkg_versions)[idx])
}

## If FRK is already present but not FRK version 2.0.1, 
## add it to the list of packages to be installed. 
if (!("FRK" %in% new_packages) && packageVersion("FRK") < '2.0.1') 
  new_packages <- c(new_packages, "FRK")


## Now install the new packages: Here, we always install the correct package 
## version (why not)
if(length(new_packages)) {
  cat("Package dependencies are being installed automatically using scripts/Dependencies_install.R\n")
  for (pkg in new_packages) {
    devtools::install_version(pkg, version = pkg_versions[pkg], repos = CRANMIRROR)
  }
}

## Change the already installed packages to the correct versions IF we have been told to do so
if(length(already_installed_pkgs_different_versions)) {
  if (exists("install_correct_versions") && install_correct_versions) {
    for (pkg in already_installed_pkgs_different_versions) {
      devtools::install_version(pkg, version = pkg_versions[pkg], repos = CRANMIRROR)
    }
  } 
} 


# if(length(new_packages)) {
#   warning("Package dependencies are being installed automatically 
#            via scripts/Dependencies_install.R")
# 
#   if (exists("install_correct_versions") && install_correct_versions) {
#     if (!("devtools" %in% rownames(installed.packages()))) {
#       install.packages("devtools", repos = CRANMIRROR)
#     }
#     for (pkg in new_packages) {
#       devtools::install_version(pkg, version = pkg_versions[pkg], repos = CRANMIRROR)
#     }
#   } else {
#     install.packages(new_packages, repos = CRANMIRROR)
#   }
#   
# } else {
#   message("\nNo package dependencies to be installed.\n")
# }