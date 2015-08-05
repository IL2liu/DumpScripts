#!/bin/sh

####################################################################
# rename files with only their id                                  #
# 3000bin-hg19-H1-204.fastq.map.bed into 204.bed                   #
####################################################################

for filename in *.bed
do
   echo $filename
   newname=`echo $filename | awk 'BEGIN {FS="-"} {print substr($4,1,3)}'` 
   newfile=$newname".bed"                       
   echo $newfile
   mv $filename $newfile
done
