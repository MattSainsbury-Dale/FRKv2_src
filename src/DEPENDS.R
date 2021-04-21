if (!("renv" %in% rownames(installed.packages())))
  install.packages("renv")
DEPENDS <- renv::dependencies()
DEPENDS <- unique(DEPENDS$Package)
new_packages <- DEPENDS[!(DEPENDS %in% rownames(installed.packages()))] 

## If FRK is already present but not FRK v.2, 
## add it to the list of packages to be installed. 
if (!("FRK" %in% new_packages) && packageVersion("FRK") < 2)
  new_packages <- c(new_packages, "FRK")

if(length(new_packages)) install.packages(new_packages)

write.table(DEPENDS, "DEPENDS.txt", row.names = FALSE, col.names = FALSE)


