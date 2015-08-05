#!/usr/bin/env python
# -*- coding: utf-8 -*-

###########################################################################
# get list of files without FAIL in                                       #
# all the fastqc quality check                                            #
# e.g ./get_FilesWithoutFail.py 'output/ALLSummary.txt' > noFail.txt     #
###########################################################################

import sys
from pprint import pprint
from collections import defaultdict

dict_grad = defaultdict(list)
all_files = []

with open(sys.argv[1]) as f_:
   line_no = 0
   for line in f_:
      items = line.split('\t')
      #criteria: PASS/FAIL/WARN
      criteria = items[0]
      file_n = items[2].strip('\n')
      all_files.append((file_n,criteria))
      for k,v in all_files:
         dict_grad[k].append(v)
      line_no += 1

for file_ in dict_grad:
   n_fail = dict_grad[file_].count("FAIL")
   # get only the files without FAIL
   if n_fail == 0:
      sys.stdout.write("%s\tNO FAILs\n"%(file_))
