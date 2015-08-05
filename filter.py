#!/usr/bin/env python
# -*- coding:utf-8 -*-

##########################################
# get the lines that are not duplicated  #
##########################################

import sys

with open(sys.argv[1]) as f:
   previousline = f.readline()
   previous = previousline.rstrip().split('\t')[0]
   for line in f:
      item = line.rstrip().split('\t')
      if item[0] != previous:
         sys.stdout.write(previousline)
         previousline = line
         previous = item[0]
      else:
         previousline = ''

sys.stdout.write(previousline)

