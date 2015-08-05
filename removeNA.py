#!/usr/bin/env python
# -*- coding: utf-8 -*-

#date created: 031012
#date modified: 051012
 #get the header (first line)

#######################################################
# remove all the lines if all the columns are NA      #
#######################################################

import sys

input_f = sys.argv[1]

with open(input_f) as f:
   for row in f:
      # get the number of reads
      reads = row.split()[3::]
      if reads.count('NA') == len(reads):
         #skip if all the reads are zero
         continue
      sys.stdout.write(row)
