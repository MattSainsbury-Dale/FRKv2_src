filenames=results/* 
for file in $filenames; do
  extension="${file#*.}"
  if [ $extension == "html" ]; 
    then cat $file >> all.html; 
  elif [ $extension == "png" ];
    # then echo "<img src=$file width=700 height=600> <br>" >> all.html
    then echo "<img src=$file> <br>" >> all.html
  fi; 
done

# Use the file number to add section titles (couldn't figure this out)
# sec=0
# filenames=`ls results`
# current_sec=$file | cut -c1-3 | sed -e s/[^0-9]//g #| awk '$1=$1' FS= OFS="."
  # if [ "$(($current_sec > $sec))" ];
  # # if [ $current_sec > $sec ];
  #   then sec=$currentsec; echo "Section title $current_sec" >> all.html;
  # fi;