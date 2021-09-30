#!/bin/bash


filenames=results/*
currentsecnum=0
for file in $filenames; do 
  
  # Do we need to add a section heading?
  name=${file##*/}                                          # remove path from file
  secnum=`echo $name | cut -c1-3  | sed -e s/[^0-9]//g`     # extract the first 3 characters, then keep only the numbers (i.e., drop the underscore _)
  if [ $secnum -gt $currentsecnum ]; then                   # Note: Cannot use > here: See https://unix.stackexchange.com/a/353633
  currentsecnum=$secnum; 
  secnum=`echo $secnum | awk '$1=$1' FS= OFS="."`           # Add "." between digits
  echo "<br> <H1> Section $secnum <br>" >> all.html;        # Write the heading
  fi;

  # Add the file depending on its extension:
  extension="${file#*.}"
  if [ $extension == "html" ]; 
    then echo "<br> <object data=$file width=700></object> <br>"  >> all.html; 
  elif [ $extension == "png" ];
    then echo "<img src=$file width=700> <br>" >> all.html
  fi; 
done