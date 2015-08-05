#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys

from d_go import d_go
from d_uniqueColor50rule import d_uniqueColor50rule

# get the unique gene associated with at least one GO id
uniqueGene = set([gene for genelist in d_go.values() for gene in genelist])

# color list
color_list = [d_uniqueColor50rule[i] if i in d_uniqueColor50rule else "NA" for i in uniqueGene]
# count the coverage of each go_id in each color
color_count = dict(zip(color_list, map(color_list.count, color_list)))
sys.stdout.write("%d\t%d\t%d\t%d\t%d\n"%( color_count['Red'],
                                          color_count['Pink'],
                                          color_count['Yellow'],
                                          color_count['Black'],
                                          color_count['NA']))
