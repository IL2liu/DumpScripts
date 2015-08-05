#/bin/sh

trap 'echo Keyboard interruption... ; exit 1' SIGINT

if [ $# -eq 0 ]; then
  echo "Usage -i [input_file] -o [output_file]"
  echo "[input_file] : list of .fastq files"                             
  echo "[output_file}: list of .fastq files!"
exit 1
fi

while getopts ":i:o:" opt; do
  case $opt in
    i)
      echo "-i your input file is: $OPTARG" >&2
      input_f="$OPTARG"
      ;;
    o)
      echo -e "-o output file: $OPTARG\n" >&2
      output_f="$OPTARG"
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

touch $output_f
for i in $(ls $input_f)
do
  num_reads=`awk '{sum +=$4} END {print sum}' $i`;
  f_name=${i##*/};
  echo $f_name$'\t'$num_reads >> $output_f
done
