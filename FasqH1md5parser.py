#!/usr/bin/python
#date of birth: 191012

#script to get the md5 information from the annotation files

import re

def H1md5Parser(annot_file, filename):
   f = open(annot_file, 'r')
   with open(filename , 'w') as fh:
      for row in f:
         if row[0] == '#':
            continue

         #get the addInfo records
         addInfo = row.split(',')[8].split('|')
         addInfo = ','.join(addInfo)

         #get md5
         md5 = re.search('md5sum=.*', addInfo).group(0)[7::]

         #get filename
         fname = str(row.split(',')[10])

         #md5 check with md5 and filename
         md5check =md5 + '  ' + fname

         fh.write(md5check)

if __name__ == '__main__':
   #get the md5 fasq-H1-hESC and stored it as FasqH1md5Source.txt
   H1md5Parser('Fasq-H1-hESC_annotation.txt', 'FasqH1md5Source.txt')
