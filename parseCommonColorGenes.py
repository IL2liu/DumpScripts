#! /usr/bin/env python
# -*- coding: utf-8 -*-

# parse the commond ids that correspond to the same color
# e.g
# NSG00000238587,ENSG00000244583 Black , parsed into
# NSG00000238587 Black
# ENSG00000244583 Black

import sys
import re

with open(sys.argv[1]) as f:
   for line in f:
      items = line.rstrip().split('\t')
      ids = items[0]
      color = items[1]
      # only if there is a comma between id, e.g
      # ENSG00000238587,ENSG00000244583
      if re.search(',', ids):
         common_ids = ids.split(',')
         for i in common_ids:
            sys.stdout.write("%s\t%s\n"%(i, color))
      else:
            sys.stdout.write("%s\t%s\n"%(ids, color))
