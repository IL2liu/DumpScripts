#!/bin/sh
# date of created: 051212

############################################################
# filter out the reads that were map with multiple matches #
############################################################

for i in $(cat $@):
  base_name=${i##*/};
  sort -k1 $i | uniq -u > "Output_dir/$base_name" && echo "OK" >&2;
done
