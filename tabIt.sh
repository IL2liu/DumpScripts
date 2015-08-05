#!/bin/bash

trap 'echo cruel interruption;exit 1' SIGINT

for i in $(ls *.bed)
do
  echo $i
  awk -v OFS="\t" '$1=$1' $i > tmp
  mv tmp $i
done

