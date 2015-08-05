#!/usr/bin/env python
# -*-: coding: utf-8 -*-

import sys

with open(sys.argv[1], 'r') as f:
   for line in f:
      items = line.rstrip().split('\t')
      if items[3] == '0:1':
         sys.stdout.write(line)
      continue


