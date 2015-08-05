#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import sys
from collections import defaultdict
import operator
import itertools
import re

with open(sys.argv[1], 'r') as in_f:
   virus_freq = defaultdict(int)
   for line in in_f:
      items = line.rstrip().split('\t')
      virus_info = items[0].split('|')
      virus_id = virus_info[1]
      virus_name =  re.search(r' ([^:]+)', virus_info[4]).groups()[0]
      virus_freq[virus_id, virus_name] += 1

   sortedVirus_freq=sorted(virus_freq.iteritems(), key = operator.itemgetter(1), reverse=True)
   total_freq = len(sortedVirus_freq)
   for (virus_id, virus_name), freq in sortedVirus_freq:
      sys.stdout.write('%s\t%s\tfreq:%d\tscore:%.5f\n'%(virus_id, virus_name, freq, freq/total_freq))

