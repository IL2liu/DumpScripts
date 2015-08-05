#!/usr/bin/env python
# -*- coding:utf-8- -*-

import sys

S = set([])
# Method 2.
# hashSet = {}
with open(sys.argv[1]) as f:
   for line in f:
      item = line.split('\t')
      # Make sure the sequence is unique.
      if item[1] in S: continue
      # Method 2.
      # if hashSet.get(item[0], False): pass
      # hashSet[item[0]] = True
      S.add(item[1])
      sys.stdout.write(line)
