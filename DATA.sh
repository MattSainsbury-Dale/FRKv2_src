if [ -e  data/chicago_crime_df.Rda ] 
then 
  echo "The Chicago crime data has already been downloaded."
else 
  curl -o data/chicago_crime_df.Rda https://dl.dropboxusercontent.com/s/cwlziqavrpe7ode/chicago_crime_df
fi 

if [ -e  data/Sydney_shapefiles/SA1/SA1_2011_AUST.shp ] 
then 
  echo "The Sydney SA1 shapefiles have already been downloaded."
else 
  curl -o data/Sydney_shapefiles/SA1/SA1_2011_AUST.shp https://dl.dropboxusercontent.com/s/v75fpk0r6pig8o4/SA1_2011_AUST.shp
fi 

if [ -e  data/Sydney_shapefiles/SA2/SA2_2011_AUST.shp ] 
then 
  echo "The Sydney SA2 shapefiles have already been downloaded."
else 
  curl -o data/Sydney_shapefiles/SA2/SA2_2011_AUST.shp https://dl.dropboxusercontent.com/s/m5ocngvauxosz50/SA2_2011_AUST.shp
fi 

if [ -e  data/Sydney_shapefiles/SA3/SA3_2011_AUST.shp ] 
then 
  echo "The Sydney SA3 shapefiles have already been downloaded."
else 
  curl -o data/Sydney_shapefiles/SA3/SA3_2011_AUST.shp https://dl.dropboxusercontent.com/s/m4qigm9s4ijn4bj/SA3_2011_AUST.shp
fi 
