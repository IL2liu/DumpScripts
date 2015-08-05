#!/usr/bin/python
# -*- coding:utf-8 -*-

import re
import sys
import tempfile
from gzopen import gzopen

BYTES = 65536

s = int(sys.argv[1])
pad = ''.join(['#']*36)

with tempfile.TemporaryFile() as temp:
   # Make a temp fasta file without newline on the sequence.
   with gzopen(sys.argv[2]) as f:
      txt = f.read(BYTES)
      while txt != '':
         while '>' in txt:
            txt += f.read(BYTES)
            header = re.search(r'\n?>[^\n]+\n', txt)
            temp.write(txt[:header.start()].replace('\n', ''))
            temp.write(txt[header.start():header.end()])
            txt = txt[header.end():]
         temp.write(txt.replace('\n', ''))
         txt = f.read(BYTES)

   # Reset temp file and read line by line.
   temp.seek(0)
   for line in temp:
      if line[0] == '>':
         seqname = line[1:].rstrip()
         continue
      for i in xrange(0, len(line)-36, s):
         sys.stdout.write('@%s_%d\n%s\n+\n%s\n' %
            (seqname, i+1, line[i:(i+36)], pad))
