#/bin/sh 
#date created: 191012

#extract all the .fastq files from the .tgz files

trap 'echo Keyboard interruption... ; exit 1' SIGINT

while getopts ":i:o:" opt; do
  case $opt in
    i)
      echo "-i your input files are: $OPTARG" >&2
      input="$OPTARG"
      ;;
    o)
      echo "-o your output file is stored: $OPTARG" >&2
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


for fastq in $(ls $input)
do                                                                     
  input_f=${fastq##*/} 
  echo "currently processing: $input_f"
  #grep the list of fastq files
  tar tvfz "$fastq" --wildcards '*/*.fastq' | awk '{print $6}' > tmp && echo "Listing File OK" >&2;
  #start extracting the list of files
  for i in $(cat tmp)
  do
    tar xzf "$fastq" "$i" -C . && echo "Extracting file OK" >&2;
    f_name="$fastq-$(basename $i)" 
    mv "$i" ./"$f_name"
  done
done                           

