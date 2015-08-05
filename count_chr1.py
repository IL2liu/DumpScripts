#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys
from collections import defaultdict

counts = defaultdict(int)

coord_max = 0
with open("chr1.txt") as f:
   for line in f:
      items = line.split()
      a = int(items[2]) / 40000
      b = int(items[5]) / 40000
      counts[(a,b)] += 1
      if a != b:
         counts[(b,a)] += 1
      if a > coord_max: coord_max = a
      if b > coord_max: coord_max = b

for i in range(coord_max+1):
   line = '\t'.join([str(counts[(i,j)]) for j in range(coord_max+1)])
   sys.stdout.write(line + "\n")
