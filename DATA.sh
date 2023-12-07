#!/bin/bash

## NB: Comments are the commands I used when I was using Dropbox to store the data sets, retained just in case I need them. 

if [ -e  data/chicago_crime_df.Rda ] 
then 
  echo "The Chicago crime data has already been downloaded."
else 
  curl -o data/chicago_crime_df.Rda https://zenodo.org/records/10278749/files/chicago_crime_df.Rda?download=1
fi 

if [ -e  data/Sydney_shapefiles/SA1 ] 
then 
  echo "The Sydney SA1 shapefiles have already been downloaded."
else 
  curl -o SA1.zip https://zenodo.org/records/10278720/files/SA1.zip?download=1
  unzip SA1.zip
  mv SA1 data/Sydney_shapefiles
  rm SA1.zip
  rm -r __MACOSX
fi 

if [ -e  data/Sydney_shapefiles/SA2 ] 
then 
  echo "The Sydney SA2 shapefiles have already been downloaded."
else 
  curl -o SA2.zip https://zenodo.org/records/10278720/files/SA2.zip?download=1
  unzip SA2.zip
  mv SA2 data/Sydney_shapefiles
  rm SA2.zip
  rm -r __MACOSX
fi 

if [ -e  data/Sydney_shapefiles/SA3 ] 
then 
  echo "The Sydney SA3 shapefiles have already been downloaded."
else 
  curl -o SA3.zip https://zenodo.org/records/10278720/files/SA3.zip?download=1
  unzip SA3.zip
  mv SA3 data/Sydney_shapefiles
  rm SA3.zip
  rm -r __MACOSX
fi 
