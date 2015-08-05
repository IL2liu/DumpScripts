#/bin/sh -i
#date of birth: 170912

#############################################
# create only filename without the content  #
# used to run test the renaming script      #
#############################################

for long_filename in $(ls Output_Data/*.bed)
do
  #remove the parent directory
  short_filename=${long_filename##*/};
  touch $short_filename 
done
