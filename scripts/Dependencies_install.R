if (!("devtools" %in% rownames(installed.packages()))) {
  install.packages("devtools", repos = CRANMIRROR)
}

CRANMIRROR <- "https://cran.csiro.au" 

tmp <- read.table("dependencies.txt", header = FALSE)
pkg_versions <- setNames(as.character(tmp[, 2]), tmp[, 1])
rm(tmp)


# ---- Non-CRAN packages ----

## These packages are treated individually because they are not available on 
## CRAN, so we need to specify their repos. 

## dggrids is simple to deal with, so we just install the current version: 
if(!("dggrids" %in% rownames(installed.packages())))
  install.packages("dggrids", repos="https://andrewzm.github.io/dggrids-repo", type= "source")


if(!("INLA" %in% rownames(installed.packages()))) {
  if (exists("install_exact_versions") && install_exact_versions) {
    # devtools::install_version("INLA", 
    #                           # repos = "https://inla.r-inla-download.org/R/", 
    #                           repos = "https://inla.r-inla-download.org/R/stable", 
    #                           version = pkg_versions["INLA"])
    
    ## Can't get the above to work. Not sure how to download exact versions from 
    ## non-standard repos. Just installing the current stable versions for now. 
    install.packages("INLA", repos="https://inla.r-inla-download.org/R/stable")
    
  } else {
    install.packages("INLA", repos="https://inla.r-inla-download.org/R/stable")
  }
}

## Remove these packages from the search list so that the script does not 
## attempt to re-install them
pkg_versions <- pkg_versions[!(names(pkg_versions) %in% c("INLA", "ddgrids"))]

# ---- CRAN packages ----

## Find the packages not yet installed and add them to the list
installed_idx <- names(pkg_versions) %in% rownames(installed.packages())
new_packages  <- names(pkg_versions)[!(installed_idx)] 

if (exists("install_exact_versions") && install_exact_versions) {
  ## Find the packages that are installed, but not the correct version
  installed_pkg_versions <- sapply(names(pkg_versions)[installed_idx], function(pkg) as.character(packageVersion(pkg)))
  idx          <- installed_pkg_versions != pkg_versions[installed_idx]
  already_installed_pkgs_different_versions <- names(installed_pkg_versions)[idx]
}

## If FRK is already present but not FRK version 2.0.1, 
## add it to the list of packages to be installed. 
if (!("FRK" %in% new_packages) && packageVersion("FRK") < '2.0.1') 
  new_packages <- c(new_packages, "FRK")

## Now install the new packages: Here, we always install the correct 
## package version (no reason not to)
if(length(new_packages)) {
  cat("Package dependencies are being installed automatically using scripts/Dependencies_install.R\n")
  for (pkg in new_packages) {
    devtools::install_version(pkg, version = pkg_versions[pkg], repos = CRANMIRROR)
  }
}

## Change the already installed packages to the correct versions IF we have been told to do so
if(exists("install_exact_versions") && install_exact_versions && length(already_installed_pkgs_different_versions)) {
    for (pkg in already_installed_pkgs_different_versions) {
      devtools::install_version(pkg, version = pkg_versions[pkg], repos = CRANMIRROR)
  } 
} 


