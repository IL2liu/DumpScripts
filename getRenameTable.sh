#/bin/sh 

#for safety
trap 'echo Keyboard interruption... ; exit 1' SIGINT

if [ $# -eq 0 ]; then
  echo "./getRenameTable.sh [input files] -o [output]"
  echo "e.g ./getRenameTable.sh -i annot.txt -o renameTable.txt"
exit 1
fi                                                                            

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

# get uncommented file
awk '{ if(substr($1, 0, 1) != "#") {print $0} }' $input > tmp1
# get only the download links and new names 
awk -F',' '{print $9"\t"$10}' tmp1 > tmp2  
# get only the download names
awk -F'/' '{print $8}' tmp2 > $output
rm tmp1 tmp2
