#!/usr/bin/env python
# -*- coding:utf-8 -*-

#date of birth: 190912

###############################################################
# Rename None assembly file names that have the known GEO no  #
###############################################################

import os

#dict for mapping the reference assembly
#mapping after manually check the GEO page on pubmed
dict_none={
      'None-H1-510.bed': 'hg19-H1-510.bed',
      'None-H1-512.bed': 'hg19-H1-512.bed',
      'None-H1-522.bed': 'hg19-H1-522.bed',
      'None-H1-524.bed': 'hg19-H1-524.bed',
      'None-H1-538.bed': 'hg19-H1-538.bed',
      'None-H1-542.bed': 'hg19-H1-542.bed',
      'None-H1-544.bed': 'hg19-H1-544.bed',
      'None-H1-546.bed': 'hg19-H1-546.bed',
      'None-H1-548.bed': 'hg19-H1-548.bed',
      'None-H1-550.bed': 'hg19-H1-550.bed',
      'None-H1-552.bed': 'hg19-H1-552.bed',
      'None-H1-554.bed': 'hg19-H1-554.bed',
      'None-H1-556.bed': 'hg19-H1-556.bed',
      'None-H1-558.bed': 'hg19-H1-558.bed',
      'None-H1-560.bed': 'hg19-H1-560.bed',
      'None-H1-564.bed': 'hg19-H1-564.bed',
      'None-H1-566.bed': 'hg19-H1-566.bed',
      'None-H1-568.bed': 'hg19-H1-568.bed',
      'None-H1-572.bed': 'hg19-H1-572.bed',
      'None-H1-574.bed': 'hg19-H1-574.bed',
      'None-H1-584.bed': 'hg19-H1-584.bed',
      'None-H1-586.bed': 'hg19-H1-586.bed',
      'None-H1-600.bed': 'hg19-H1-600.bed',
      'None-H1-602.bed': 'hg19-H1-602.bed',
      'None-H1-604.bed': 'hg19-H1-604.bed',
      'None-H1-606.bed': 'hg19-H1-606.bed',
      'None-H1-608.bed': 'hg19-H1-608.bed',
      'None-H1-610.bed': 'hg19-H1-610.bed',
      'None-H1-612.bed': 'hg19-H1-612.bed',
      'None-H1-614.bed': 'hg19-H1-614.bed',
      'None-H1-616.bed': 'hg19-H1-616.bed',
      'None-H1-620.bed': 'hg19-H1-620.bed',
      'None-H1-622.bed': 'hg19-H1-622.bed',
      'None-H1-624.bed': 'hg19-H1-624.bed',
      }

#start renaming ...
for filename in os.listdir("Output_Data_3000bpBin/."):
   if filename in dict_none.keys():
      new_name = dict_none[filename]
      os.rename(os.path.join("Output_Data_3000bpBin", filename),
            os.path.join("Output_Data_3000bpBin", new_name))
      continue
