#/bin/sh

###############################################################################################################
# Checked the simple stats of raw_fastq files,                                                                #
# script provided by [GGD]                                                                                    #
# (http://gettinggeneticsdone.blogspot.com.es/2012/04/awk-command-to-count-total-unique-and.html),many thanks #
# e.g output:                                                                                                 #
#     99115 60567 61.1078 ACCTCAGGA 354 0.357161                                                              #
#     This is telling you:                                                                                    #
#     The total number of reads (99,115).                                                                     #
#     The number of unique reads (60,567).                                                                    #
#     The frequency of unique reads as a proportion of the total (61%).                                       #
#     The most abundant sequence (useful for finding adapters, linkers, etc).                                 #
#     The number of times that sequence is present (354).                                                     #
#     The frequency of that sequence as a proportion of the total number of reads (0.35%).                    #
###############################################################################################################

trap 'echo Keyboard interruption... ; exit 1' SIGINT

if [ $# -eq 0 ]; then
  echo "getSimpleStats.sh -i [input_dir] -o [output_file]"
  echo "e.g ./getSimpleStats.sh -i input -o 'output/ALLSummary.txt'"
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

for read in `ls $input/*.gz`; 
  do echo -n "$read "; 
  gunzip -c $read | awk '((NR-2)%4==0){read=$1;total++;
  count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};
  if(count[read]==1){unique++}};
  print total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}'; 
  done > $output
