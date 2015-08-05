#!/bin/sh

trap 'echo Keyboard interruption... ; exit 1' SIGINT

array=( chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr1 chr20 chr21 chr22 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chrX chrY )

for chr in "${array[@]}"
do
  paste *_$chr.txt > tmp.txt
  echo -e "H1\t$chr" > header.txt
  cat header.txt tmp.txt > "H1_"$chr"_binary".txt
  rm header.txt
done
