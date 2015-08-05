####################################################
# Pipeline to get the intersect of two bed files   #
####################################################

trap 'echo Keyboard interruption... ; exit 1' SIGINT

if [ $# -eq 0 ]; then
  echo "Usage -a [input_fileA] -b [list_fileB] -s [suffixName] -o [output_dir]"
  echo "[input_fileA] : input_file to be intersected"
  echo "[list_fileB]: features to be intersected [list of bed files stored in !one directory!]"
  echo "[suffixName]: suffix name to be inserted in the output, eg .h1hESCColor.bed"
  echo "[output_dir]: storage directory !create this dir in advance!"
  echo "./getBedIntersect.sh -a test/NoHeaderStates_dom_hESC.txt.bed -b test_COORDS/ -s .test -o test/"
exit 1                                                                                        
fi

while getopts ":a:b:s:o:" opt; do
  case $opt in
    a)
      echo "-a your input fileA is: $OPTARG" >&2
      input_fA="$OPTARG"
      ;;
    b)
      echo "-b features to be intersected in: $OPTARG" >&2
      input_fB="$OPTARG"
      ;;
    s)
      echo "-s suffix to be inserted: $OPTARG" >&2
      suffixName="$OPTARG"
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

for i in $(ls $input_fB)
do
  f_=`echo $i | awk -F'.' '{print $1}'`
  output_name=$f_"$suffixName"
  intersectBed -a $input_fA -b $input_fB/$i > $output_dir/$output_name 
done


