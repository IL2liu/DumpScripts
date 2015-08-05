#/bin/sh

##################################
# combined all the summary files #
##################################

trap 'echo Keyboard interruption... ; exit 1' SIGINT

if [ $# -eq 0 ]; then
  echo "Usage:getEncodingFastqC.sh -i [input_dir] -o [output_file]"
  echo "e.g ./getEncodingFastqC.sh -i output/ -o 'output/ALLSummary.txt'"
exit 1
fi

while getopts ":i:o:" opt; do
  case $opt in
    i)
      echo "-i your input dir is: $OPTARG" >&2
      input="$OPTARG"
      ;;
    o)
      echo "-o your output file: $OPTARG" >&2
      output="$OPTARG"
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
list_dir=`ls -l $input | egrep '^d' | awk '{print $9}'`

touch $output
for c_dir in $list_dir
do
  #get the encoding
  encode=$(sed -n 6p $input/$c_dir/fastqc_data.txt)
  echo -e "$c_dir\t$encode" >> $output
done

