#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import sys
from collections import defaultdict
import operator
import itertools

##########################
# Filter the virus reads #     
##########################

with open(sys.argv[1], 'r') as in_f:
   virus_freq = defaultdict(int)
   for line in in_f:
      items = line.rstrip().split('\t')
      if items[4] != '-':
         viruses = items[4].split(',')
         for virus in viruses:
            virus_freq[virus] += 1


   #sortedVirus_freq=sorted(virus_freq.iteritems(), key = operator.itemgetter(1), reverse=True)
   for item in sorted(virus_freq.iteritems(), key = operator.itemgetter(1), reverse=True):
      sys.stdout.write('%s\t%d\n'%(item))

