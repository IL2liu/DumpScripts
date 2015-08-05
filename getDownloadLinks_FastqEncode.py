#!/usr/bin/python

#script to get the download links from the annotation file

import re
import sys

def getDownloadLinks(annot_file, filename):
   f = open(annot_file, 'r')
   with open(filename, 'w') as fh:
      for row in f:
         if row[0] == '#':
            continue

         #get the download link
         link = str(row.split(',')[8]) + '\n'
         fh.write(link)

if __name__ == '__main__':
   getDownloadLinks(sys.argv[1], sys.argv[2])
