#! /usr/bin/env python
# -*- coding: utf-8 -*-

####################################################
# note that the geneID was associated with TSS     #
# each gene could be assigned with several TSS,    #
# therefore it's associated with different colors  #
####################################################

from collections import defaultdict
import sys
import operator

frequencies = defaultdict(int)
dict_gene = defaultdict(list)

with open(sys.argv[1]) as f:
   for line in f:
      items = line.rstrip().split(' ')
      gene_id = items[0]
      color = items[1]
      dict_gene[gene_id].append(color)

# counting
for ids in dict_gene:
   #import pdb; pdb.set_trace()
   colors = dict_gene[ids]
   colors
   # map the gene_id to its respective color
   # count the coverage of each go_id in each color
   color_count = dict(zip(colors, map(colors.count, colors)))
   # get the majority color
   major_color = max(color_count.iteritems(), key=operator.itemgetter(1))[0]
   sys.stdout.write("%s\t%s\n"%(ids, major_color))
