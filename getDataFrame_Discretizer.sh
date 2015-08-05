#/bin/sh
# date created: 131212

#################################################
# Combine sample with input for HMM discretizer #
#################################################

if [ $# -eq 0 ]; then                                                        
  echo "Usage:./getDataFrame_Discretizer.sh -i [input_file] -c [input_sample] -o [output_dir]"
  echo "./getDataFrame_Discretizer.sh -i 'test/*.bed' -c 'test/input.txt' -o foo &"
  echo "!mkdir output_dir FIRST!"
exit 1
fi

while getopts ":i:c:o:" opt; do
  case $opt in
    i)
      echo "-i your input files are: $OPTARG" >&2
      input="$OPTARG"
      ;;
    c)
      echo "-c your input samples is $OPTARG" >&2
      input_sample="$OPTARG"
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

for i in $(ls $input)
do
  file_n=$(basename $i)
  # get the file number
  f_no=`echo $file_n | sed 's/.*H1-\([0-9]*\).*/\1/'`

  echo -e "chrName\tstart\tend\t$f_no\tinput" > $output_dir/headerfile.txt 
  
  # combine files
  paste $i $input_sample | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' > $output_dir/tmp.txt

  # add header
  cat $output_dir/headerfile.txt $output_dir/tmp.txt > $output_dir/tmpfile.txt 

  output_f="combined-"$file_n
  mv $output_dir/tmpfile.txt $output_dir/$output_f
  rm $output_dir/tmp.txt
  rm $output_dir/headerfile.txt
done

