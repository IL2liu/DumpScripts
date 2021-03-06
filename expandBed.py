#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

with open(sys.argv[1]) as f:
   for line in f:
      item = line.rstrip().split()
      # each line looks like (below)
      # chr1  100110000   100119000   E1
      expand_range = range(int(item[1]), int(item[2]), 3000)
      for i in expand_range:
         sys.stdout.write( "%s\t%d\t%d\t%s\n"%(item[0],i+1, i+3000, item[3]) )
