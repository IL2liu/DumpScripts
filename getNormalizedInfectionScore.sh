#! /bin/sh                                                                               

########################################################################
# Pipeline to combine the binary files from ChromHMM into one big file #
########################################################################

trap 'echo Keyboard interruption... ; exit 1' SIGINT

if [ $# -eq 0 ]; then
  echo "Usage:./getNormalizedInfectionScore.sh -i [input_file] -o [output_dir]"
  echo "e.g ./getNormalizedInfectionScore.sh -i TEST -o ."
  echo "input_file : list of .map files containing reads that map to viral genomes"
  echo "output_dir: the name of your output dir to store"
exit 1
fi

while getopts ":i:o:" opt; do
  case $opt in
    i)
      echo "-i your input files are: $OPTARG" >&2
      input_f="$OPTARG"
      ;;
    o)
      echo "-o your output file is stored: $OPTARG" >&2
      output_dir="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

for i in $(cat $input_f)
do 
  # filter the map files to get unique maps 
  f_file="filtered-"$i
  ./filterVirusesRead.py $i > $f_file 
  # count the number of infections
  c_file="countedInfection-"$i
  ./countInfections.py $f_file > $c_file
  # normalized the score
  n_file="NormalizedInfection-"$i
  ./getNormalizedScore.py $c_file > $n_file
  sort -r -n -k2 $n_file > "sortedNormalizedInfection-"$i
done


