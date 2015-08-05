#! /usr/bin/env python
# -*- coding: utf-8 -*-

#####################################################
# map the gene (ENSEMBL_ID) to its respective color #
#####################################################

import sys

# dictionary containing the mapping
# d_color ={'ENSG00000227999': 'Black', ...}
from d_color import d_color

with open(sys.argv[1]) as f:
   # skip the first line as a header
   f.readline()
   for line in f:
      items = line.rstrip().split('\t')
      genes = items[5].split(', ')
      color_list = [d_color[i] if i in d_color else "NA" for i in genes]
      color_count = dict(zip(color_list, map(color_list.count, color_list)))
      sys.stdout.write("%s\t%d\t%d\t%d\t%d\n"%(items[1], color_count['Red'],
                                          color_count['Pink'], 
                                          color_count['Yellow'],
                                          color_count['Black']))
      #sys.stdout.write ("%s\t%s\n" % (items[1], ','.join(colors)))
