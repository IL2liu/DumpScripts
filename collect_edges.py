#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import sys
from collections import defaultdict
import operator
import itertools
import re

#################################################
# collect edges, pairs of possible viruses      #
#################################################

with open(sys.argv[1], 'r') as in_f:
   collect_edges= defaultdict(int)

   for line in in_f:
      items = line.rstrip().split('\t')
      # What Guillaume says:
      # sort_virus_name = sorted(re.findall(',? ([^:]+)', items[4]))
      if items[4] != '-':
         # get the list of viruses
         viruses = items[4].split(',')
         # extract only the viral names
         virus_name=[re.search(r' ([^:]+)', virus).groups()[0] for virus in viruses]
         sort_virus_name = sorted(virus_name)

         #get all possible pairs
         virus_pairs = itertools.combinations(sort_virus_name,2)
         for virus_pair in virus_pairs:
            collect_edges[virus_pair] += 1
   sorted_collect_edges = sorted(collect_edges.iteritems(), key = operator.itemgetter(1), reverse=True)

   for pair, freq in sorted_collect_edges:
      sys.stdout.write('%s\t%d\n'%(pair,freq))
   #for item in sorted(collect_edges.items()):
   #   sys.stdout.write('%s\t%d\n'% item)




