#/bin/sh                                                                               
                                                                                       
#for safety                                                                            
trap 'echo Keyboard interruption... ; exit 1' SIGINT                                   
                                                                                       
if [ $# -eq 0 ]; then                                                                  
  echo "Usage:gzipThreads.sh -i [input files] -o [output_dir] -n [threads]"      
  echo "e.g ./gzipThreads.sh -i './*.fastq' -o . -n 2"                          
exit 1                                                                                 
fi                                                                                     
                                                                                       
while getopts ":i:o:n:" opt; do                                                        
  case $opt in                                                                         
    i)                                                                                 
      echo "-i your input files are: $OPTARG" >&2                                      
      input="$OPTARG"                                                                  
      ;;                                                                               
    o)                                                                                 
      echo "-o your output directory: $OPTARG" >&2                                     
      output="$OPTARG"                                                                 
      ;;                                                                               
    n)                                                                                 
      echo "-n no.threads: $OPTARG" >&2                                                
      thread="$OPTARG"                                                                 
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
                                                                                       
count=0                                                                                
for f_ in $(ls $input)                                                              
do                                                                                     
  gzip -9 $f_ &
  let count+=1; [[ $((count%$thread)) -eq 0 ]] && wait                                 
done                                                                                   
