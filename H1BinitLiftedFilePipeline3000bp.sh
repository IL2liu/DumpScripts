#/bin/sh -i
#date of birth : 051012 

################################################################
# Pipeline to bin only lifted bed files  with the size of 3000 #
################################################################

trap 'echo Keyboard interruption... ; exit 1' SIGINT

# Run through the bed files of the directory.
# Output_Data_300bpBin directory contains the previosly lifted bed file
for long_bedfile in $(ls Output_Data_300bpBin/*.bed)
do
  short_bedfile=${long_bedfile##*/};
  if [[ "${short_bedfile:0:6}" == "lifted" ]]; then
    echo "binning $short_bedfile..." >&2;

    #the binned file is name without lifted, thus as follow: hg19-H1-***.bed
    binned_bed=${short_bedfile/lifted-/}

    #start binning here (bin size is 300)
    ./binit.py -b 3000 -f bed $long_bedfile > "Output_Data_3000bpBin/$binned_bed" && echo "Went OK" >&2;
  fi  
done
