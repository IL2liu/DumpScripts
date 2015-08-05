#! /usr/bin/env python
# -*- coding: utf-8 -*-

# cols 1-3 are the genomic coordinates: chr start end
# cols 4. AT content
# cols 5. low complexity seq

import sys

with open(sys.argv[1], "r") as f:
   f.readline()
   seq = f.read()

seq = seq.replace("\n", '')
for i in xrange(len(seq)/3000):
   chunck = seq[i*3000:(i+1)*3000]
   N = chunck.count('N')
   AT = chunck.count('A') + chunck.count('T')
   at = chunck.count('a') + chunck.count('t')
   gc = chunck.count('g') + chunck.count('c')
   sys.stdout.write('chr1\t%d\t%d\t' % (1+i*3000, (i+1)*3000))
   try:
      sys.stdout.write('%f\t%f\n' % \
            (float(at+AT)/ (3000-N), float(at+gc)/ (3000-N)))
   except ZeroDivisionError:
      sys.stdout.write('NA\tNA\n')
