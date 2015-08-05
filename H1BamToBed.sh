#/bin/sh -i
#date of birth : 270912 

#######################################################################
# Pipeline to convert the hg18 bam files to hg18 bed file             #
# for the end position add arbitrarily 100 nucleotides from the start #
#######################################################################

trap 'echo Keyboard interruption... ; exit 1' SIGINT

# Run through the bam files of the directory.
for long_bamfile in $(ls Input_Data/*.bam)
do
  short_bam=${long_bamfile##*/}; 
  if [[ "${short_bam:0:4}" == "hg18" ]]; then
    echo "converting ${short_bam} from hg18.bam to hg18.bed ...." >&2; 
    #change the file extension from bam to bed
    hg18_bed="${short_bam%bam}bed"
    samtools view $long_bamfile | awk -v OFS="\t" '{print $3, $4, $4+100}' > "Output_Data/$hg18_bed" && echo "OK" >&2;
  fi
done
