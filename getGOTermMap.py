#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
import re
from vtrack import vskip

##########################
# get the mapping table  #
# GOID   GOTERM          #

# get the dict_go from G.F
dict_go = json.load(vskip(open(sys.argv[1])))

for i in dict_go.keys():
   goID = re.match('^GO:[0-9]+',i).group(0)
   goTERM = re.search('\((.*)\)', i).group(1)
   sys.stdout.write("%s\t%s\n"%(goID, goTERM))
