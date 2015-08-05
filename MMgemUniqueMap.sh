#/bin/sh 

###############################################
# map fastq files                             #
# !for mouse genome!                          #
# gem version: 1.366                          # 
###############################################

trap 'echo Keyboard interruption... ; exit 1' SIGINT

if [ $# -eq 0 ]; then
  echo "Usage -i [input_list] -n [no.core] -m [mismatches] -o [output_dir]"
  echo "[input_list] : list of .fastq files"
  echo "[no.core]: for threading"
  echo "[mismatches]: no.of allowed mismatches"
  echo "[output_dir}: storage directory !create this dir in advance!"
exit 1
fi

while getopts ":i:n:m:o:" opt; do
  case $opt in
    i)
      echo "-i your input file is: $OPTARG" >&2
      input_f="$OPTARG"
      ;;
    n)
      echo "-n number of cores uses: $OPTARG" >&2
      core="$OPTARG"
      ;;
    m)
      echo "-n number of mismatches allowed: $OPTARG" >&2
      mismatch="$OPTARG"
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
  filename="${i_file%.*.*}"
  gunzip -c $i | gem-mapper -I '/software/private/gfilion/seq/gem/mm10/mm10_masked.gem' -q 'ignore' -m $mismatch --unique-mapping -o $output_dir/$filename -T $core && echo "Gem mapping of $i:OK" >&1;
done                                            
