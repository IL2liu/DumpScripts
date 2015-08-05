#!/usr/bin/env python
# -*- coding: utf-8 -*-

# date : 031112

###############################################################
# create the chromHMM cellmarkfiletable from annotation file  #
###############################################################

import re
import argparse
import sys

def cellmarkfiletable(input_annotFile, output_name):

   with open(output_name , 'w') as output_f:
      for line in input_annotFile:
         if line[0] == "#":
            continue
         else:
            items = line.split(',')
            add_infos = items[8].split('|')
            cell = add_infos[0]
            cell_name = re.sub('cell=','', cell)
            mark = add_infos[2]
            mark_name = re.sub('antibody=','',mark)
            file_name = "%s.bed" % (items[10].strip('\n'))
            output_f.write("%s\t%s\t%s\n" %  \
                  (cell_name,mark_name,file_name))

if __name__ == "__main__":
   # Build argument parser.
   parser = argparse.ArgumentParser(
         prog = 'cellmarkfiletable.py',
         description = """create the cellmarkfiletable to run chromHMM"""
   )
   parser.add_argument(
      'input_file',
      metavar = 'f',
      type = file,
      nargs = '?',
      default = sys.stdin,
      help = 'file to process'
   )
   parser.add_argument(
      '-o',
      '--output_name',
      metavar = 'filename',
      type= str,
      nargs='?',
      help = 'name for the outputfile'
   )
   args = parser.parse_args()

   cellmarkfiletable(args.input_file, args.output_name)
