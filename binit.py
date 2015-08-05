#!/usr/bin/env python
# -*- coding:utf-8 -*-
import sys
import os
from collections import defaultdict

def binning(dir_bam):
   """
   sync the bin index for all the bam files
   """

   path = '../%s'%(dir_bam)
   for bamfile in glob.glob(os.path.join(path, '*.bam')):
      bamfile = bamfile.split('/')[-1]
      print bamfile
      bam_txt = bamfile+".txt"
      os.system('samtools view %s > %s'%(bamfile, bam_txt)
      
   
   """
   hist = defaultdict(int)

   for line in sys.stdin:
      items = line.split()
      hist[(items[2],int(items[3])/300)] += 1

   chroms = set([chrom for (chrom,bin) in hist.keys()])

   for chrom in chroms:
      maxbin = max([bin for (chrom__,bin) in hist.keys() if chrom__ == chrom])
      for bin in range(maxbin):
         sys.stdout.write("%s\t%s\t%s\t%d\n" % \
               (chrom, 1+bin*300, (1+bin)*300, hist[(chrom,bin)]))
   """

if __name__=="__main__":
   binnit('cells')
