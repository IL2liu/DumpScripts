#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import gzip
from gzopen import gzopen

with open(sys.argv[1]) as in_f:
   l = []
   for line in in_f:
      l.append(line)

out_f = gzip.open(sys.argv[2], 'w')
for i in l:
   out_f.write(i)
out_f.close()


