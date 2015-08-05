#/bin/sh

for i in $(ls *.bed)
do
  sum_=`awk '{ sum+=$4} END {print sum}' $i`
  if [[ $sum_ < 0 ]]
  then 
    echo $i
  fi
done
