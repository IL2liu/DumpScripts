#/bin/sh -i
#date of birth: 190912

#############################################
# Scipt to check the lifted over sequences  #               
#                                           #
#############################################

trap 'echo Keyboard interruption... ; exit 1' SIGINT

for long_liftedBed in $(ls Output_Data/*.bed)
do
  short_liftedBed=${long_liftedBed##*/};
  if [[ "${short_liftedBed:0:11}" == "lifted-hg19" ]]; then  
    echo "Check lifting over of $short_liftedBed...." >&2

    hg18Bed=${short_liftedBed/lifted-hg19/hg18}     #get the original file on hg18 that was used for lifting over
    unliftedCo_hg18="unliftedComment-"$hg18Bed      #get the unlifted file containing comments           
    unlifted_hg18="unlifted-"$hg18Bed;

    #grep only the coordinates from unmmaped file 
    sed -e '/^#/d' "Output_Data/$unliftedCo_hg18" > "Output_Data/$unlifted_hg18"

    # remove the coordinate of hg18 that were not lifted -over due to deletions
    lifted_hg18="lifted-"$hg18Bed
    diff "Output_Data/$hg18Bed" "Output_Data/$unlifted_hg18" | grep "<" | sed -e 's/^< //' > "Output_Data/$lifted_hg18" 
    
    #check with wc-l that the total number of lines between hg18-H1-002Lifted.bed should be the same with hg19-H1-002.bed
    lineNo_hg18=(wc -l "Output_Data/$lifted_hg18") | awk '{print $1;}'
    lineNo_hg19=(wc -l "Output_Data/$short_liftedBed") | awk '{print $1;}'
 
    if [[ $lineNo_hg18 == $lineNo_hg19 ]]; then
      echo "Number of lines before and after lifting over is equal...sync went well"  >&2
      #run the script to randomize the bed coordinates for checking
      #note adjust the dir_store on the pickRandomBedCoordinate script where to store ... must be in the same as in the bash script
      #e.g in this case the dir_store is in Output_Data
      ./pickRandomBedCoordinate.py "Output_Data/$short_liftedBed" "Output_Data/$lifted_hg18"
      #change the format for the coordinates in order to use twoBitToFa script
      pickRandomBedAfter="randomlyPicked-"$short_liftedBed
      pickRandomBedBefore="randomlyPicked-"$lifted_hg18                                                 

      #change the format of bed file for compatibility with twoBitToFa script
      coordinateAfter="coordinate-"$pickRandomBedAfter
      coordinateBefore="coordinate-"$pickRandomBedBefore
      awk '{print $1":"$2"-"$3}' "Output_Data/$pickRandomBedAfter"  > "Output_Data/$coordinateAfter"   
      awk '{print $1":"$2"-"$3}' "Output_Data/$pickRandomBedBefore"  > "Output_Data/$coordinateBefore"

      #change .bed into .fa
      FastaPickRandomBedAfter=${pickRandomBedAfter%bed}fa
      FastaPickRandomBedBefore=${pickRandomBedBefore%bed}fa
      twoBitToFa -seqList="Output_Data/$coordinateAfter" hg19.2bit "CheckLiftOver/$FastaPickRandomBedAfter"                                         
      twoBitToFa -seqList="Output_Data/$coordinateBefore" hg18.2bit "CheckLiftOver/$FastaPickRandomBedBefore"

      diff --ignore-case "CheckLiftOver/$FastaPickRandomBedAfter" "CheckLiftOver/$FastaPickRandomBedBefore"  > "CheckLiftOver/diff-${short_liftedBed%bed}txt" &&  echo "checking went OK" >&2;

    else
      echo "Number of lines are not sync, please check the files $hg18Bed" >&2
    fi

  fi   
done

