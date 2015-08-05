#!/bin/bash

if [ $# -eq 0 ]; then                                                        
  echo "Usage:./chromHMMparse.sh -i [input_file] -o [output_dir]"
  echo "./chromHMMparse.sh -i 'test/*.bed' -o foo &"
  echo "!mkdir output_dir FIRST!"
exit 1
fi

while getopts ":i:c:o:" opt; do
  case $opt in
    i)
      echo "-i your input files are: $OPTARG" >&2
      input="$OPTARG"
      ;;
    o)
      echo "-od your output file is stored: $OPTARG" >&2
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


trap 'echo Keyboard interruption... ; exit 1' SIGINT

for f in $(ls $input)
do
  echo $f
  b_name=${f%.bed}
  id=${b_name##*/}
  for chr in `cut -f 1 $f | sort | uniq`; 
  do
    output_f=$output_dir/$id"_"$chr.txt
    grep -w $chr "$f" | awk -F'\t' '{print $4}' > $output_f
    echo $id > header.txt
    cat header.txt $output_f > tmp
    mv tmp $output_f
  done
done
