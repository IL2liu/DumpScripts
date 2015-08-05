#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys

from d_go import d_go
from d_uniqueColorMajority import d_uniqueColorMajority

# get the unique gene associated with at least one GO id
uniqueGene = set([gene for genelist in d_go.values() for gene in genelist])

# color list
color_list = [d_uniqueColorMajority[i] if i in d_uniqueColorMajority else "NA" for i in uniqueGene]
# count the coverage of each go_id in each color
color_count = dict(zip(color_list, map(color_list.count, color_list)))
sys.stdout.write("%d\t%d\t%d\t%d\t\n"%( color_count['Red'],
                                          color_count['Pink'],
                                          color_count['Yellow'],
                                          color_count['Black'],
                                          ))
                                          
