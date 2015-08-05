# -*- coding:utf-8 -*-

import os
import re
import sys
from itertools import izip

from gzopen import gzopen

def merge(fnamelist):
   flist = [gzopen(fname) for fname in fnamelist]

   # Get names and print header.
   idexpr = r'-(\d{3}[a-z]?)'
   IDs = [re.search(idexpr, fname).group(1) for fname in fnamelist]
   sys.stdout.write('seqname\tstart\tend\t' + '\t'.join(IDs) + '\n')

   # Iterate through all the files at the same time with 'izip'.
   for linetuple in izip(*flist):
      # Extract seqname, start and end from first file.
      mapping = '\t'.join(linetuple[0].split()[:3])
      # Extract 4-th column and print.
      entries = '\t'.join([line.split()[3] for line in linetuple])
      sys.stdout.write(mapping + '\t' + entries + '\n')

   for f in flist:
      f.close()


if __name__ == '__main__':
   merge(sys.argv[1:])
