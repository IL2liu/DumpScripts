#! /usr/bin/env python
# -*- coding: utf-8 -*-

#####################################################
# get the binary matrix of GO and genes association #
# matrix looks like this
#   genes...(cols)
# G
# O
# .
# .
# (GO IDs in the rows)
#####################################################

import sys
# d_go dictionary from amigo
from d_go import d_go

# get all the unique genes
all_genes = set([gene for genes in d_go.values() for gene in genes])

# print the header
# Genes_ID ..... GO_ID
sys.stdout.write('\t'.join(all_genes) + '\tGO_ID\n')
for go_id in d_go:
   for gene in all_genes:
      sys.stdout.write('%d\t'%(0+ (gene in d_go[go_id])))
   sys.stdout.write('%s\n'%(go_id))
