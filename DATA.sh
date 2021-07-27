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

if [ -e  data/Sydney_shapefiles/SA1 ] 
then 
  echo "The Sydney SA1 shapefiles have already been downloaded."
else 
  #curl -o data/Sydney_shapefiles/SA1/SA1_2011_AUST.shp https://dl.dropboxusercontent.com/s/v75fpk0r6pig8o4/SA1_2011_AUST.shp
  wget -np -nH --cut-dirs 5 https://hpc.niasra.uow.edu.au/ckan/dataset/9559f1ba-bd8b-41a7-8229-726e7bf88db8/resource/8300cace-6885-47fd-a021-37595ac8dc84/download/sa1.zip
  unzip sa1.zip
  rm sa1.zip
  mv SA1 data/Sydney_shapefiles
fi 

if [ -e  data/Sydney_shapefiles/SA2 ] 
then 
  echo "The Sydney SA2 shapefiles have already been downloaded."
else 
  #curl -o data/Sydney_shapefiles/SA2/SA2_2011_AUST.shp https://dl.dropboxusercontent.com/s/m5ocngvauxosz50/SA2_2011_AUST.shp
  wget -np -nH --cut-dirs 5 https://hpc.niasra.uow.edu.au/ckan/dataset/9559f1ba-bd8b-41a7-8229-726e7bf88db8/resource/5390f101-8636-472c-b1dc-f7066e9f9898/download/sa2.zip
  unzip sa2.zip
  rm sa2.zip
  mv SA2 data/Sydney_shapefiles
fi 

if [ -e  data/Sydney_shapefiles/SA3 ] 
then 
  echo "The Sydney SA3 shapefiles have already been downloaded."
else 
  # curl -o data/Sydney_shapefiles/SA3/SA3_2011_AUST.shp https://dl.dropboxusercontent.com/s/m4qigm9s4ijn4bj/SA3_2011_AUST.shp
  wget -np -nH --cut-dirs 5 https://hpc.niasra.uow.edu.au/ckan/dataset/9559f1ba-bd8b-41a7-8229-726e7bf88db8/resource/59607912-e088-4429-a096-90638c1159b0/download/sa3.zip
  unzip sa3.zip
  rm sa3.zip
  mv SA3 data/Sydney_shapefiles
fi 
