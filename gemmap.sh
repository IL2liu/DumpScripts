#/bin/sh 

###############################################
# map fastq files                             #
# gem version: 1.366                          # 
###############################################

INDEX=/software/private/gfilion/seq/gem/hg19/hg19_unmasked.gem
QUALITY=ignore
MM=2
THREADS=4

queue=.gem_queue.ls
output_dir=.


while getopts "i:o:" opt; do
  case $opt in
    i)
      i="$OPTARG"
      ;;
    o)
      output_dir="$OPTARG"
      ;;

    \?)
      echo "invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "option -$OPTARG requires an argument." >&2
      exit 1                                                
      ;;
  esac
done


if [ ! -f "$queue" ]; then
  echo "creating queue from file $i" >&2
  cat $i >> $queue
elif [ "X$i" != "X" ]; then
  echo "** queue exists: ignoring input file **" >&2
fi
     
# 'mkdir' is atomic, so we can use it for locking.
lock=/tmp/gemmap_queue.lock

trap 'echo -e "Keyboard interruption \n** $geminput **" ; exit 1' SIGINT
trap 'rm -rf "$lock"' EXIT

while [ -f "$queue" ]; do
  if mkdir "$lock"; then
    geminput=$(head -1 $queue)
    sed -i '1d' $queue
    count=$(wc -l $queue)
    if [ "$count" == "0 $queue" ]; then
      rm -rf $queue
      break
    fi
    rm -rf "$lock"
    echo "processing $geminput" >&2
    gemoutput=${geminput##.*/}
    gemoutput=${gemoutput%.gz}
    gemoutput=${gemoutput%.fastq}
    gunzip -c $geminput | gem-mapper -I $INDEX -q $QUALITY -m $MM \
      --unique-mapping -o $output_dir/$gemoutput -T $THREADS 2>/dev/null
  else
    sleep 1
  fi
done

exit 0


#if [ $# -eq 0 ]; then
#  echo "Usage -i [input_list] -n [no.core] -m [mismatches] -o [output_dir]"
#  echo "[input_list] : list of .fastq files"
#  echo "[no.core]: for threading"
#  echo "[mismatches]: no.of allowed mismatches"
#  echo "[output_dir]: storage directory !create this dir in advance!"
#exit 1
#fi
#
#while getopts ":i:n:m:o:" opt; do
#  case $opt in
#    i)
#      echo "-i your input file is: $OPTARG" >&2
#      input_f="$OPTARG"
#      ;;
#    n)
#      echo "-n number of cores used: $OPTARG" >&2
#      core="$OPTARG"
#      ;;
#    m)
#      echo "-m number of mismatches allowed: $OPTARG" >&2
#      mismatch="$OPTARG"
#      ;;
#    o)
#      echo -e "-o output directory: $OPTARG\n" >&2
#      output_dir="$OPTARG"
#      ;;
#    \?)
#      echo "Invalid option: -$OPTARG" >&2
#      exit 1
#      ;;
#    :)
#      echo "Option -$OPTARG requires an argument." >&2
#      exit 1                                                
#      ;;
#  esac
#done
# 
