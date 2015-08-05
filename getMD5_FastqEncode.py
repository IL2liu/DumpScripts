#!/usr/bin/python

#script to get the md5 information from the annotation files

import sys
import re

def md5Parser(annot_file, filename):
   f = open(annot_file, 'r')
   with open(filename , 'w') as fh:
      for row in f:
         if row[0] == '#':
            continue

         #get the addInfo records
         addInfo = row.split(',')[7].split('|')
         addInfo = ','.join(addInfo)

         #get md5
         md5 = re.search('md5sum=.*', addInfo).group(0)[7::]

         #get filename
         fname = str(row.split(',')[8].split('/')[-1])

         #md5 check with md5 and filename
         md5check =md5 + '  ' + fname + '\n'

         fh.write(md5check)

if __name__ == '__main__':
   md5Parser(sys.argv[1], sys.argv[2])
