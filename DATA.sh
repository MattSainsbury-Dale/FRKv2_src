## NB: Comments are the commands I used when I was using Dropbox to store the data sets, retained just in case I need them. 

if [ -e  data/chicago_crime_df.Rda ] 
then 
  echo "The Chicago crime data has already been downloaded."
else 
  #curl -o data/chicago_crime_df.Rda https://dl.dropboxusercontent.com/s/cwlziqavrpe7ode/chicago_crime_df
  wget -np -nH --cut-dirs 5 https://hpc.niasra.uow.edu.au/ckan/dataset/9aa14ca4-796a-41a6-b1f0-03758ed5ec9a/resource/f7d9ea8f-1506-4046-bf98-a2424935a3bc/download/chicagocrimedf.rda
  mv chicagocrimedf.rda chicago_crime_df.Rda        # rename for consistency
  mv chicago_crime_df.Rda data/chicago_crime_df.Rda # move to /data
fi 

if [ -e  data/Sydney_shapefiles/SA1/SA1_2011_AUST.shp ] 
then 
  echo "The Sydney SA1 shapefiles have already been downloaded."
else 
  #curl -o data/Sydney_shapefiles/SA1/SA1_2011_AUST.shp https://dl.dropboxusercontent.com/s/v75fpk0r6pig8o4/SA1_2011_AUST.shp
  wget -np -nH --cut-dirs 5 https://hpc.niasra.uow.edu.au/ckan/dataset/9559f1ba-bd8b-41a7-8229-726e7bf88db8/resource/e021e308-f3e3-4f3b-b4c0-32809ffd72c5/download/sa12011aust.zip
  unzip sa12011aust.zip 
  rm sa12011aust.zip 
  mv SA1_2011_AUST.shp data/Sydney_shapefiles/SA1/SA1_2011_AUST.shp 
fi 

if [ -e  data/Sydney_shapefiles/SA2/SA2_2011_AUST.shp ] 
then 
  echo "The Sydney SA2 shapefiles have already been downloaded."
else 
  #curl -o data/Sydney_shapefiles/SA2/SA2_2011_AUST.shp https://dl.dropboxusercontent.com/s/m5ocngvauxosz50/SA2_2011_AUST.shp
  wget -np -nH --cut-dirs 5 https://hpc.niasra.uow.edu.au/ckan/dataset/9559f1ba-bd8b-41a7-8229-726e7bf88db8/resource/af61a4e6-8ced-4364-9f43-937cd15fceea/download/sa22011aust.zip
  unzip sa22011aust.zip 
  rm sa22011aust.zip 
  mv SA2_2011_AUST.shp data/Sydney_shapefiles/SA2/SA2_2011_AUST.shp 
fi 

if [ -e  data/Sydney_shapefiles/SA3/SA3_2011_AUST.shp ] 
then 
  echo "The Sydney SA3 shapefiles have already been downloaded."
else 
  # curl -o data/Sydney_shapefiles/SA3/SA3_2011_AUST.shp https://dl.dropboxusercontent.com/s/m4qigm9s4ijn4bj/SA3_2011_AUST.shp
  
  wget -np -nH --cut-dirs 5 https://hpc.niasra.uow.edu.au/ckan/dataset/9559f1ba-bd8b-41a7-8229-726e7bf88db8/resource/13795155-ee32-428c-ba04-2417ee210ab7/download/sa32011aust.zip
  unzip sa32011aust.zip 
  rm sa32011aust.zip 
  mv SA3_2011_AUST.shp data/Sydney_shapefiles/SA3/SA3_2011_AUST.shp 
fi 
