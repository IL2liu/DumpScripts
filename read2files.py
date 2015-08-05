#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys
from gzopen import gzopen

def testquality(s, t):
   for char in s:
      if ord(char) < 63: return False
   for char in t:
      if ord(char) < 63: return False
   return True

def test(s, t):
   # Return 'True' if 's' and 't' are close to 'P1' and 'P2',
   # resp., 'False' otherwise.
   P1 = 'GTATTTGGTAGCATTGCCTTT'
   P2 = 'AAACTTCCGACTTCAACTGT'
   i = 0
   j = 0
   for k in xrange(20):
      if P1[k] != s[k]:
         i += 1
         if i > 3: return False
      if P2[k] != t[k]:
         j += 1
         if j > 3: return False
   if (i == 3 and P1[20] != s[20]): return False
   return True


def match_reads(file_name_1, file_name_2):
   outf = {}
   passed = False
   file1 = gzopen(file_name_1)
   file2 = gzopen(file_name_2)
   for (lineno, line1) in enumerate(file1):
      line2 = file2.readline()
      modulo = lineno % 4
      if modulo == 0 and passed:
         try:
            outf[index].write('@' + spotname + seq + '+\n' + quality)
         except KeyError:
            outf[index] = open(index + '.fastq', 'w')
            outf[index].write('@' + spotname + seq + '+\n' + quality)
      elif modulo == 1:
         index = line1[3:7]
         passed = test(line1[7:], line2)
         if passed:
            seq = line2[20:]
            spotname = line1[28:]
      elif modulo == 3 and passed:
         passed = testquality(line1[28:50], line2[20:50])
         quality = line2[20:]

   file1.close()
   file2.close()
   for key in outf: outf[key].close()

if __name__ == '__main__':
   match_reads(sys.argv[1], sys.argv[2])
