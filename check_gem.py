#!/usr/bin/python
# -*- coding:utf-8 -*-

import sys
from gzopen import gzopen

OK = 0
NO = 0
total = 0
with gzopen(sys.argv[1]) as f:
   for line in f:
      items = line.split()
      from_ = items[0].split('_')
      to_ = items[4].split(':')
      if from_[0] == to_[0] and from_[1] == to_[2]:
         OK += 1
      elif items[4] != '-':
         NO += 1
      total += 1

sys.stdout.write('total correct %d (%.2f%%)\n'
      'total incorrect %d (%.2f%%)\n' %
      (OK, 100*float(OK)/total, NO, 100*float(NO)/total))
