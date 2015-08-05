#!/usr/bin/python
# -*- coding:utf-8 -*-

import sys
from gzopen import gzopen

OK = 0
NO = 0
total = 0
with gzopen(sys.argv[1]) as f:
   # Skip header (and sacrifice first line).
   for line in f:
      if line[0] != '@': break
   for line in f:
      items = line.split()
      from_ = items[0].split('_')
      to_ = (items[2], items[3])
      if from_[0] == to_[0] and from_[1] == to_[1]:
         OK += 1
      elif items[1] != '*':
         NO += 1
      total += 1

sys.stdout.write('total correct %d (%.2f%%)\ntotal incorrect %d (%.2f)\n' %
      (OK, 100*float(OK)/total, NO, 100*float(NO)/total))
