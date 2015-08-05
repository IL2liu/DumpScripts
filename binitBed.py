#!/usr/bin/env python
# -*- coding:utf-8 -*-


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

# default assembly: hg19 reference genome.
limits_hg19 = {
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

limits_mm10 = {
    'chr1': 195471971,
    'chr2': 182113224,
    'chr3': 160039680,
    'chr4': 156508116,
    'chr5': 151834684,
    'chr6': 149736546,
    'chr7': 145441459,
    'chr8': 129401213,
    'chr9': 124595110,
    'chr10': 130694993,
    'chr11': 122082543,
    'chr12': 120129022,
    'chr13': 120421639,
    'chr14': 124902244,
    'chr15': 104043685,
    'chr16': 98207768,
    'chr17': 94987271,
    'chr18': 90702639,
    'chr19': 61431566,
    'chrX': 171031299,
    'chrY': 91744698,
}

def binit(bin_size, input_file, limits=limits_hg19, output_dir='.'):
   """
   bin into certain window size
   """

   getmap = itemgetter(0,2)
   unique_maps = set()
   unique_counts = defaultdict(int)
   multiple_counts = defaultdict(int)

   with gzopen(input_file) as f:
      for line in f:
         # Fields: read name, sequence, quality, map count, positions(s).
         item = line.rstrip().split('\t')
         # keep only reads with unique map
         # the following scenarios are accepted
         # 1.. ; 0:0:...:0:1...; 0:0:0+1... 
         if item[4] == '-': continue
         stratum_size = int(re.search(r'^[0:+]*(\d+)', item[3]).groups()[0])
         #match_uniq = re.search(r'^[0:+]*1(?:\Z|\D)', item[3])
         thisdict = unique_counts if stratum_size == 1 else multiple_counts
         for hit in item[4].split(',')[:stratum_size]:
            #mapping = item[4]
            # 'mapping' is like "chr1:+:12942:34T1,chr15:-:102518193:34T1"
            #chrom,pos = getmap(mapping.split(':'))
            chrom,pos = getmap(hit.split(':'))
            # Keep only one read (the first) of a read series.
            if (chrom,pos) in unique_maps: continue
            unique_maps.add((chrom,pos))
            thisdict[(chrom,int(pos)/bin_size)] += 1

   #output file
   if not os.path.exists(output_dir):
      try:
         os.makedirs(output_dir)
      except OSError as exception:
         if exception.errno != errno.EEXIST:
            raise

   head, tail = os.path.split(input_file)
   base = os.path.splitext(tail)[0]
   output_fname="%sbin-%s.bed" %(bin_size,base)
   output_file = str(os.path.join(output_dir, output_fname))

   with open(output_file, 'w') as output_f:
      for chrom,size in sorted(limits_hg19.items()):
         for b in range(size/bin_size):
            coord = (chrom,b)
            mapping = "%s\t%d\t%d\t" % (chrom, 1+b*bin_size, (1+b)*bin_size)
            unmappable = multiple_counts[coord] and not unique_counts[coord]
            count = 'NA' if unmappable else str(unique_counts[(chrom,b)])
            output_f.write(mapping + count + '\n')


class binThread (multiprocessing.Process):
   """
   Binning in threads
   """
   def __init__(self, lock, queue, bin_size, output_dir):
      self.lock = lock
      self.queue = queue
      self.bin_size = bin_size
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

         binit(bin_size=self.bin_size,
               output_dir=self.output_dir,\
               input_file=input_f)
         self.queue.task_done()


if __name__ == "__main__":
   # Build argument parser.
   parser = argparse.ArgumentParser(
         prog = 'binitGemMapThread.py',
         description = """ Binning the aligment file in 300bp(default bin size):
         !use hg19 reference genome!
         - n: number of cores for multiprocessing
         !NOTE this program requires a big chunk of memory, 
         the binning dumped directly into memory to speed up the binning!""",
         formatter_class=RawTextHelpFormatter
   )
   parser.add_argument(
         '-r',
         '--ref-limit',
         metavar = 'r',
         type = str,
         nargs = '?',
         default = 'limits_hg19',
         help = 'chromosomes sizes of the reference genome'
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
   threads = [binThread(lock, queue, args.bin_size, args.ref_limit, args.output_dir) for i in range(args.n_threads)]
   for thread in threads:
      thread.start()
   queue.join()

   # Fill in task queue while threads are running.
   for f_task in glob.glob(args.input_file):
      queue.put_nowait(f_task)
