#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
      sort_virus_name = sorted(re.findall(',? ([^:]+)', items[4]))

      #get all possible pairs
      virus_pairs = itertools.combinations(sort_virus_name,2)
      for virus_pair in virus_pairs:
         collect_edges[virus_pair] += 1

   # sort based on the frequency 
   #for item in sorted(collect_edges.iteritems(), key = operator.itemgetter(1), reverse=True):
   for item in sorted(collect_edges.items()):
      sys.stdout.write('%s u %s\n'% (item[0][0].replace(" ","_"), item[0][1].replace(" ", "_")))
