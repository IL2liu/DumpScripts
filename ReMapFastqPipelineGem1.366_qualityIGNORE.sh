#/bin/sh 

###############################################
# Pipeline to convert .fastq to .sam with GEM #
# gem version: 1.366                          # 
###############################################

trap 'echo Keyboard interruption... ; exit 1' SIGINT

if [ $# -eq 0 ]; then
  echo "Usage -i [input_file] -n [no.core] -o [output_dir]"
  echo "[input_file] : list of .fastq files"
  echo "[no.core]: for threading"
  echo "[output_dir}: storage directory !create this dir in advance!"
exit 1
fi

while getopts ":i:n:o:" opt; do
  case $opt in
    i)
      echo "-i your input file is: $OPTARG" >&2
      input_f="$OPTARG"
      ;;
    n)
      echo "-n number of cores uses: $OPTARG" >&2
      core="$OPTARG"
      ;;
    o)
      echo -e "-o output directory: $OPTARG\n" >&2
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
  #mapping ...
  i_file=${i##*/}; 
  #echo "Gem-Mapping:$i_file"
  gem-mapper -I '/software/private/gfilion/seq/gem/hg19/hg19_masked.gem' -i $i -q 'ignore' -o $output_dir/$i_file -T $core && echo "Gem mapping of $i:OK" >&1;
  map_f=$i_file.map
done                                            
