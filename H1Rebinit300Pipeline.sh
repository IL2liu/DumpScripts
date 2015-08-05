#/bin/sh -i
#date of birth : 170917 

###################################################################
# Pipeline to rebin only the .bam files that were not OK          #
# List of error files was stored in Checks/H1Binit300logNOTOK.txt #
###################################################################

trap 'echo Keyboard interruption... ; exit 1' SIGINT

# Run through the bam files of the directory that are not ok
for longbamfile in $(cat Checks/H1Binit300logNOTOK.txt)
do
  echo "processing ${longbamfile%bam}bed..." >&2;
  bamfile=${longbamfile##*/};
  samtools view $longbamfile | ./binit.py -b 300 > "Output_Data/${bamfile%bam}bed" && echo "OK" >&2;
done
