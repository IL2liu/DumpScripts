#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
with open(sys.argv[1]) as f:
   for line in f:
      if line[0:3] =="chr":
         items = line.split(":")
         sys.stdout.write("%s:%s\n"%(items[0], items[2]))
      else:
         sys.stdout.write(line)
