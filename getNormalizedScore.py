#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

# import dict containing the name of virus(key) and size(value)
from dict_virusIDSize import dict_virusIDSize


import sys

with open(sys.argv[1], "r") as in_f:
   for line in in_f:
      items = line.rstrip().split('\t')
      v_id = items[0]
      v_name = items[1].replace(" ","_")
      v_freq = int(items[2])
      v_size = int(dict_virusIDSize[v_id])
      v_Nfreq = v_freq/v_size
      sys.stdout.write("%s\t%f\n" %(v_name, v_Nfreq))
