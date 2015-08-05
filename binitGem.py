#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys
import gc
import re

from operator import itemgetter
from collections import defaultdict
from gzopen import gzopen

# default assembly: hg19 reference genome.
limits_hg19 = {
   'chr1': 249250621,
   'chr2': 243199373,
   'chr3': 198022430,
   'chr4': 191154276,
   'chr5': 180915260,
   'chr6': 171115067,
   'chr7': 159138663,
   'chr8': 146364022,
   'chr9': 141213431,
   'chr10': 135534747,
   'chr11': 135006516,
   'chr12': 133851895,
   'chr13': 115169878,
   'chr14': 107349540,
   'chr15': 102531392,
   'chr16': 90354753,
   'chr17': 81195210,
   'chr18': 78077248,
   'chr19': 59128983,
   'chr20': 63025520,
   'chr21': 48129895,
   'chr22': 51304566,
   'chrX': 155270560,
   'chrY': 59373566,
}

def binit(bin_size, fname, limits=limits_hg19):

   getmap = itemgetter(0,2)
   unique_reads = set([])
   counts = defaultdict(int)

   with gzopen(fname) as f:
      for line in f:
         # Fields: read name, sequence, quality, map count, positions(s).
         item = line.rstrip().split('\t')
         # Keep only reads with unique map: "0:1" or "1+...".
         match_uniq = re.search(r'^[0:+]*1(?:\Z|\D)', item[3])
         if match_uniq:
            mapping = item[4]
            # 'mapping' is like "chr6:+:52132829:2C29".
            chrom,pos = getmap(mapping.split(':'))
            # Keep only on read (the first) of a read series.
            if (chrom,pos) in unique_reads: continue
            unique_reads.add((chrom,pos))
            counts[(chrom,int(pos)/bin_size)] += 1
         continue


   for chrom,size in sorted(limits_hg19.items()):
      for b in range(size/bin_size):
         sys.stdout.write("%s\t%d\t%d\t%d\n" % \
            (chrom, 1+b*bin_size, (1+b)*bin_size, counts[(chrom,b)]))

if __name__ == "__main__":
   binit(3000, sys.argv[1])

