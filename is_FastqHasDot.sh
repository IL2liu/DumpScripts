#/bin/sh

##################################################################
# get the list of fastq files that have the reads with \. (dot)  #
# yes the file contains reads with '.' otherwise no              #
##################################################################

trap 'echo keyboard interruptionn ; exit 1' SIGINT

if [ $# -eq 0 ]; then
  echo "Usage: is_FastqHasDot.sh -i [input_dir]"
  echo "e.g ./is_FastqHasDot.sh -i TEST" 
exit 1
fi

while getopts ":i:" opt; do
  case $opt in
    i)
      echo "-i your input directory: $OPTARG" >&2
      input="$OPTARG"
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

for i in $(ls $input/*.gz)
do
  #get the reads only
  gunzip -c $i | awk '{if (NR % 4 == 2) {print $0}}' > tmp.txt
  # test if it contains \.
  if grep -q '\.' tmp.txt ; then echo $i":yes"; else echo $i":no"; fi
  rm tmp.txt
done
