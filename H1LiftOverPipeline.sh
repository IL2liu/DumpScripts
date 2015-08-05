#/bin/sh -i
#date of birth : 170912 
#1st adjustment : 210912
  #Adjust the script to store the unLifted files

###########################################################
# Pipeline to liftover bed files  after binning in 300 bp #
# liftingover hg18 to hg19 assembly                       #
# output sent to Output_Liftover                          #
###########################################################

trap 'echo Keyboard interruption... ; exit 1' SIGINT

# Run through the bam files of the directory.
for long_bedfile in $(ls Output_Data/*.bed)
do
  short_bed=${long_bedfile##*/}; 
  if [[ "${short_bed:0:4}" == "lifted" ]]; then
    echo "Lifting over ${short_bed} from hg18 to hg19...." >&2; 
    #remove the parent directory of the filename
    lifted_bed="lifted-"${short_bed/hg18/hg19}
    #create the filename to store the unmapped files (the one that is not lifted over)
    unliftedCo_bed="unliftedComment-"$short_bed
    liftOver $long_bedfile hg18ToHg19.over.chain "Output_Data/$lifted_bed" "Output_Data/$unliftedCo_bed" && echo "Went OK" >&2;
  fi
done
