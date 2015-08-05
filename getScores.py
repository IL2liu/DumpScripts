#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from collections import defaultdict

# import dict containing the name of virus(key) and size(value)
from dict_virusIDSize import dict_virusIDSize
import re
import sys

with open(sys.argv[1], "r") as in_f:
   virus_score = defaultdict(float)
   for line in in_f:
      items = line.rstrip().split('\t')
      v_info = items[4]
      if v_info != '-':
         viruses = items[4].split(',')
         n_viruses = len(viruses)
         for virus in viruses:
            virus_info = virus.split('|')
            virus_id = virus_info[1]
            virus_name =  re.search(r' ([^:]+)', virus_info[4]).groups()[0]
            virus_name = virus_name.replace(' ','_')
            virus_score[virus_id, virus_name] += 1/(n_viruses)
   for item in virus_score.items():
      sys.stdout.write('%s\t%f\n'% (item[0][1], item[1]/dict_virusIDSize[item[0][0]] ))


