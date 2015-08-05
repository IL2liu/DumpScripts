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

   for item in sorted(virus_freq.iteritems(), key = operator.itemgetter(1), reverse=True):
      sys.stdout.write('%s\t%s\t%d\n'% (item[0][0], item[0][1], item[1]))

