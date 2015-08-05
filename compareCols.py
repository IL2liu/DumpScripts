#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
with open(sys.argv[1]) as f:
   for line in f:
      files = line.rstrip().split('\t')
      if (str(files[0]) != str(files[1])):
         sys.stdout.write(line)
