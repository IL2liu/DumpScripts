#!/usr/bin/env python
# -*- coding:utf-8 -*-

# date of birth : 130912
# date of modified : 280912 
#  - modified to add new feature that can bin bed file
# date of modified: 231012
#  - modified the input_file argument to allow the list of files as input
#  - modified to allow multiple threading

#######################################################
# Bin the alignment file, indexed always start from 1 #
# Bin index is the alignment coordinate/bin size      #
#######################################################

import os
import sys
#default dict use for counting
from collections import defaultdict
import argparse
from argparse import RawTextHelpFormatter
import threading
import Queue

#key:chr and value: size
dict_genome_size_hg19 = {
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

def dict_GenomeSize(input_file):
   """
   create a dictionary
   - key: chromosomes in the genome
   - value: size 
   """
   with open(input_file) as f:
      # Assume no header.
      # Assume 2 items per line!!
      return dict([line.split() for line in f])

def binit(bin_size, input_file, format='sam', dict_hg19 = dict_genome_size_hg19, output_dir='.'):
   """binning the coordinate of the alignment"""

   #key: chromosome, and bin index
   #value: read coverage
   freq_table = defaultdict(int)

   if format == 'bed':
      chrom_field = 0
      start_field = 1
   else:
      # Default.
      chrom_field = 2
      start_field = 3

   try:
      for line in input_file:
         items = line.split()
         if items[start_field] == '*':
            continue
         freq_table[(items[chrom_field],int(items[start_field])/bin_size)] += 1

   except IndexError:
      sys.stderr.write("Check if your file is properly formatted, field separator is any space.")
   except:
      sys.stderr.write('file error:%s'%(input_file.name))
      raise


   #collect all the aligned chromosome without duplicates 
   chroms = set([chrom for (chrom,bin) in freq_table.keys()])

   #output file name
   output_fname= '%sbin-%s' %(bin_size, input_file.name)
   if not os.path.exists(output_dir):
      os.makedirs(output_dir)

   output_file = str(os.path.join(output_dir, output_fname))

   #start binning ...
   with open (output_file, 'w') as output_f:
      for chrom in sorted(dict_hg19):
         chrom_size = int(dict_hg19[chrom])
         maxbin = int(chrom_size/bin_size)
         for bin in range(maxbin):
            output_f.write("%s\t%s\t%s\t%d\n" % \
                        (chrom, 1+bin*bin_size,
                              (1+bin)*bin_size, freq_table[(chrom,bin)]))

         """
            sys.stdout.write("%s\t%s\t%s\t%d\n" % \
               (chrom, 1+bin*bin_size,
               (1+bin)*bin_size, freq_table[(chrom,bin)]))
         """


# threading this job
class myThread (threading.Thread):
   def __init__(self, bin_size, output_dir):
      self.bin_size = bin_size
      self.output_dir = output_dir
      threading.Thread.__init__(self)

   def run(self):
      input_f = None
      while True:
         # Wait for data in task queue.
         queuelock.acquire()
         try:
            input_f = taskqueue.get_nowait()
         except Queue.Empty:
            # finish
            return
         finally:
            queuelock.release()

         # Process data.
         if input_f is not None:
            binit(bin_size=self.bin_size, output_dir=self.output_dir, input_file=input_f)
            taskqueue.task_done()
            input_f = None


         

if __name__ == "__main__":
   # Build argument parser.
   parser = argparse.ArgumentParser(
         prog = 'binit.py',
         description = """ Binning the aligment file in 300bp(default bin size), formatted as in:
         - default format is sam
           - sam files obtained by reading bam file with samtools
         - optional format is in bed """,
         formatter_class=RawTextHelpFormatter
   )
   parser.add_argument(
         '-b',
         '--bin-size',
         metavar = 'n',
         type = int,
         nargs = '?',
         default = 300,
         help = 'bin size in base pairs'
   )
   parser.add_argument(
         '-f',
         '--format',
         metavar = 'f',
         type = str,
         nargs = '?',
         default = 'sam',
         help = 'format of the input file'
   )
   parser.add_argument(
         '-n',
         '--n-threads',
         metavar = 'n',
         type = int,
         nargs = '?',
         default = 1,
         help = 'number of threads to run in parallel'
   )
   parser.add_argument(
         '-od',
         '--output_dir',
         metavar = 'f',
         type = str,
         nargs = '?',
         default = '.',
         help = 'Directory for the output file'
   )
   parser.add_argument(
         'input_file',
         metavar = 'f',
         type = file,
         nargs = '*',
         default = sys.stdin,
         help = 'file to process'
   )

   args = parser.parse_args()
 
   #start binning ....
   #binit(args.bin_size, arddgs.format, args.input_file)
 
   queuelock = threading.Lock()
   taskqueue = Queue.Queue(-1) # No limit.

   # Fill in task queue while threads are running.
   for f_task in args.input_file:
      taskqueue.put_nowait(f_task)

   # Create and start three threads.
   threads = [myThread(args.bin_size, args.output_dir) for num in range(args.n_threads)]
   for thread in threads:
      thread.start()

   # Wait for all 'task_done()' signals.
   taskqueue.join()
