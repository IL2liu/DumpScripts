#! /usr/bin/env python
# -*- coding: utf-8 -*-

######################################################
# map the gene (ENSEMBL_ID) to its respective color  #
# and count the coverage of each GO_id in each color #
# output (tab-delimited file):                       #
# go_id  Red   Pink  Yellow   Black   NA             #
######################################################

import sys

# parsed dictionary from human_go_GF20130505.json
# d_go = d_go = {"GO:0022607": ["ENSG00000102054", "ENSG00000172426",.],.} 
from d_go import d_go

# dictionary containing the mapping
# each geneID is mapped to a color (unique mapping)
# unique mapping based on the majority rule,
# as each geneID may contain several TSS, in turn, it may contain several colors
# d_uniqueColor ={'ENSG00000227999': 'Black', ...}
from d_uniqueColor import d_uniqueColor

for go_id, ensembl_id in d_go.items():
   # map the gene_id to its respective color
   color_list = [d_uniqueColor[i] if i in d_uniqueColor else "NA" for i in ensembl_id]
   # count the coverage of each go_id in each color
   color_count = dict(zip(color_list, map(color_list.count, color_list)))
   list_color = set(['NA', 'Black', 'Red', 'Yellow', 'Pink'])
   list_dict = set(color_count.keys())
   # get the color not in the list
   unmatch = list_color.symmetric_difference(list_dict)
   # put the count to zero if the color is not there
   for i in unmatch:
     color_count[i] = 0
   # import pdb; pdb.set_trace()
   sys.stdout.write("%s\t%d\t%d\t%d\t%d\t%d\n"%(go_id, color_count['Red'],
                                             color_count['Pink'], 
                                             color_count['Yellow'],
                                             color_count['Black'],
                                             color_count['NA']))
