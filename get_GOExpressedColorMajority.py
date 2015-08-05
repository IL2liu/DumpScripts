#! /usr/bin/env python
# -*- coding: utf-8 -*-

######################################################
# map the gene (ENSEMBL_ID) to its respective color  #
# and count the coverage of each GO_id in each color #
# output (tab-delimited file):                       #
# go_id  Red   Pink  Yellow   Black                  #
######################################################

import sys

# parsed dictionary from human_go_GF20130505.json
# d_go = d_go = {"GO:0022607": ["ENSG00000102054", "ENSG00000172426",.],.} 
from d_go import d_go

# dictionary containing the mapping
# d_uniqueExpressedColorMajority ={'ENSG00000227999': 'Black', ...}
from d_uniqueExpressedColorMajority import d_uniqueExpressedColorMajority

for go_id in d_go.iterkeys():
   # map the gene_id to its respective color
   color_list = [d_uniqueExpressedColorMajority[i] if i in d_uniqueExpressedColorMajority else "NA" for i in d_go[go_id]]
   # count the coverage of each go_id in each color
   color_count = dict(zip(color_list, map(color_list.count, color_list)))
   sys.stdout.write("%s\t%d\t%d\t%d\t%d\n"%(go_id, color_count['Red'],
                                             color_count['Pink'], 
                                             color_count['Yellow'],
                                             color_count['Black']))
