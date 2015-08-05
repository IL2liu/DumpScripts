#!/usr/bin/env python
# -*- coding:utf-8 -*-

####################################################################################################
# binit map file from gem-mapping                                                                  #
# gem-mapping: unique mapping with m allowed mismatches                                            #
# binit into a certain size with the size limits of each chromosome  provided from limits.py file  #
####################################################################################################

#for I/O and commandline parsing
import os
import sys
import argparse
from argparse import RawTextHelpFormatter
import glob
import gc
import errno

from operator import itemgetter
from collections import defaultdict
from gzopen import gzopen

# for threading jobs
import multiprocessing
import Queue
import limits

def binit(ref_limits, bin_size, mismatch, input_file, output_dir='.'):
   """
   select the reads that are mapped uniquely, allowing three mismatches
   bin into certain window size

   """

   try:
      # import dict containing chromosomes size
      # limits = {'chr1': 1898309, 'chr2': 2902930, ...}
      limits =__import__("limits") #limits.py is the file containing limits dict 
      limits = limits.__dict__.get(ref_limits) #choose which species, e.g hg19 , mm10
      freq_table = defaultdict(int)
   except:
      raise

   try:
      f = gzopen(input_file)
      for line in f:
         # Fields: read name, sequence, quality, map count, positions(s).
         item = line.rstrip().split('\t')
         #mapping = item[4]
         if item[4] == '-': continue
         #check for mismatches
         # e.g. 0:0:1 = allowing 2 mismatches (2 zeros)
         no_mismatch = item[3].count("0")
         if no_mismatch <= mismatch:
            # 'mapping' is like "chr1:+:12942:34T1,chr15:-:102518193:34T1"
            getmap = item[4].split(":")
            chrom = getmap[0]
            start = getmap[2]
            freq_table[(chrom, int(start)/bin_size)] += 1
   except IndexError:
      sys.stderr.write("Check if your file is properly formatted, !field separator is a tab!")
   except:
      sys.stderr.write('file error:%s'%(in_f.name))
      raise
   finally:
      f.close()

   #collect all the aligned chromosome without duplicates 
   chroms = set([chrom for (chrom,bin) in freq_table.keys()])

   #start binning ...

   # remove the garbage collector during list append
   gc.disable()

   bin_list=[]
   for chrom in sorted(limits):
      chrom_size = int(limits[chrom])
      maxbin = int(chrom_size/bin_size)
      for bin in range(maxbin):
         bin_list.append("%s\t%s\t%s\t%d\n" % \
                        (chrom, 1+bin*bin_size,
                        (1+bin)*bin_size, freq_table[(chrom,bin)]))

   gc.enable()

   #finally output the file                                     
   if not os.path.exists(output_dir):
      try:
         os.makedirs(output_dir)
      except OSError as exception:
         raise

   head, tail = os.path.split(input_file)
   base = os.path.splitext(tail)[0]
   output_fname="%sbin-%s.bed" %(bin_size,base)
   output_file = str(os.path.join(output_dir, output_fname))

   with open(output_file, 'w') as output_f:
      for line in bin_list:
         output_f.write(line)

class binThread (multiprocessing.Process):
   """
   Binning in threads
   """
   def __init__(self, lock, queue, ref_limits,  bin_size, mismatch, output_dir):
      self.lock = lock
      self.queue = queue
      self.ref_limits = ref_limits
      self.bin_size = bin_size
      self.mismatch = mismatch
      self.output_dir = output_dir
      multiprocessing.Process.__init__(self)

   def run(self):
      while True:
         # Wait for data in task queue.
         self.lock.acquire()
         try:
            input_f = queue.get_nowait()
         except Queue.Empty:
            # finish
            return
         finally:
            self.lock.release()

         binit(ref_limits = self.ref_limits,
               bin_size=self.bin_size,
               mismatch=self.mismatch,
               output_dir=self.output_dir,\
               input_file=input_f)
         self.queue.task_done()


if __name__ == "__main__":
   # Build argument parser.
   parser = argparse.ArgumentParser(
         prog = 'binitGemMap.py',
         description = """ Binning the aligment file in 300bp(default bin size):
         !use hg19 reference genome by default as a size limit for each chromosome!
         !gem unique mapping with allowed two mismatches!
         - n: number of cores for multiprocessing
         !NOTE this program requires a big chunk of memory,
         the binning dumped directly into memory to speed up the binning!""",
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
         '-m',
         '--mismatch',
         metavar = 'm',
         type = int,
         nargs = '?',
         default = 2,
         help = 'allowed # mismatches'
   )
   parser.add_argument(
         '-l',
         '--ref-limits',
         metavar = 'l',
         type = str,
         nargs = '?',
         default = 'hg19',
         help = 'dict containing chromosomes sizes of a ref.genome'
   )
   parser.add_argument(
         '-n',
         '--n-threads',
         metavar = 't',
         type = int,
         nargs = '?',
         default = 1,
         help = 'number of threads to run in parallel'
   )
   parser.add_argument(
         '-od',
         '--output_dir',
         metavar = 'od',
         type = str,
         nargs = '?',
         default = '.',
         help = 'directory for the output file'
   )
   parser.add_argument(
         'input_file',
         metavar = 'f',
         type = str,
         nargs = '?',
         default = sys.stdin,
         help = 'file to process'
   )

   args = parser.parse_args()

   # bin with threads
   lock = multiprocessing.Lock()
   queue = multiprocessing.JoinableQueue(-1) #no.limit for queue
   manager = multiprocessing.Manager()

   # Fill in task queue while threads are running.
   for f_task in glob.glob(args.input_file):
      queue.put_nowait(f_task)

   # Create and start binning in threads.

   threads = [binThread(lock, queue, args.ref_limits,  args.bin_size, args.mismatch,  args.output_dir) for i in range(args.n_threads)]
   for thread in threads:
      thread.start()
   queue.join()

   # Fill in task queue while threads are running.
   for f_task in glob.glob(args.input_file):
      queue.put_nowait(f_task)
