#/bin/sh

trap 'echo Keyboard interruption... ; exit 1' SIGINT

for i in $(ls *.bed)
do
  echo $i
  awk '{if($4 == 1) $4=0;print;}' $i | awk '{if($4 == 2) $4=1;print;}' | awk -v OFS="\t" '$1=$1' > tmp
  mv tmp $i
done
