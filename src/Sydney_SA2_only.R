library("FRK")
library("plyr")     # round_any()
library("dplyr")    
library("ggplot2")
library("sp")
library("maptools") # readShapePoly()
library("raster") # bind()
library("ggpubr")

defaultW <- getOption("warn") 
options(warn = -1)

fitting <- "SA2s"

source("./src/Sydney.R")

options(warn = defaultW)

