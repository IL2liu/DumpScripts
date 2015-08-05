#/bin/sh

#for safety
trap 'echo Keyboard interruption... ; exit 1' SIGINT

if [ $# -eq 0 ]; then
  echo "Usage:replaceDot_withN.sh -i [input list] -o [output_dir]"
  echo "e.g ./replaceDot_withN.sh -i TEST.ls -o TEST"
exit 1
fi

while getopts ":i:o:" opt; do
  case $opt in
    i)
      echo "-i your input list: $OPTARG" >&2
      input="$OPTARG"
      ;;
    o)
      echo "-o your output directory: $OPTARG" >&2
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

cat $input | while read fastq 
do
  f_=$(basename $fastq)
  f_="${f_%.*}"
  echo $f_
  gunzip -c $fastq | awk '{if (NR % 4 == 2) { gsub(/\./, "N"); print } else {print $0}}' > tmp
  mv tmp $output/$f_
done
