#/bin/sh -i
#date of birth : 140912 

#########################################################
# Pipeline to bin only .bam files  with the size of 300 #
#########################################################

trap 'echo Keyboard interruption... ; exit 1' SIGINT

# Run through the bam files of the directory.
for longbamfile in $(ls Input_Data/*.bam)
do
  echo "processing ${longbamfile%bam}bed..." >&2;
  bamfile=${longbamfile##*/};
  samtools view $longbamfile | ./binit.py -b 300 > "Output_Data/${bamfile%bam}bed" && echo "OK" >&2;
done
