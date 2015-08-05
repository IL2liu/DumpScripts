#/bin/sh

#################################################################################
# rename and create hardlink of downloaded files from ENCODE based on the list  #
#################################################################################

#for safety
trap 'echo Keyboard interruption... ; exit 1' SIGINT

if [ $# -eq 0 ]; then
  echo "Usage:renameLinkedFiles.sh -l [list_file] -i [input_dir] -o [output_dir]"
  echo "e.g ./renameLinkedFiles.sh -l list_file -i test1/ -o test2/"
exit 1
fi

while getopts ":l:i::o:" opt; do
  case $opt in
    l)
    echo "-l your list file is: $OPTARG" >&2
      list="$OPTARG"
      ;;
    i)
      echo "-i your input dir is: $OPTARG" >&2
      input="$OPTARG"
      ;;
    o)
      echo "-o your output directory is stored: $OPTARG" >&2
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

#name1:downloaded name
#name2:id name

cat $list | while read name1 name2; 
do   
  ln $input/$name1 $output/$name2; 
done
