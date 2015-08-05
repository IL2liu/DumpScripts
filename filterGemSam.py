#!/usr/bin/env python
# -*- coding:utf-8 -*-

#################################################################
# get the lines of the gemsam's output that are not duplicated  #
#################################################################

import sys

#default dict use for counting
from collections import defaultdict

#for I/O and commandline parsing
import os
import sys
import argparse
from argparse import RawTextHelpFormatter
import glob
import gc
import errno

# for threading jobs
import multiprocessing
import Queue


def filterGemSam(input_file, output_dir='.'):
 
   # remove the garbage collector during list append
   gc.disable()

   uniq_map = []
   with open(input_file) as f:
      previousline = f.readline()
      # get the map's id
      previous = previousline.rstrip().split('\t')[0]
      for line in f:
         item = line.rstrip().split('\t')
         # only if read map is not duplicated, we put into the list
         if item[0] != previous:
            uniq_map.append(previousline)
            previousline = line
            previous = item[0]
         else:
            previousline = ''
   uniq_map.append(previousline)

   gc.enable()

   #for the output the file                                     
   if not os.path.exists(output_dir):
      try:
         os.makedirs(output_dir)
      except OSError as exception:
         if exception.errno != errno.EEXIST:
            raise
   head, tail = os.path.split(input_file)
   base = os.path.splitext(tail)[0]
   output_fname="filtered-%s.sam" %(base)
   output_file = str(os.path.join(output_dir, output_fname))

   #write to the output
   with open(output_file, 'w') as output_f:
      for line in uniq_map:
         output_f.write(line)

class filterThread (multiprocessing.Process):
   """
   Filter in threads
   """
   def __init__(self, lock, queue, output_dir):
      self.lock = lock
      self.queue = queue
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

         filterGemSam(output_dir=self.output_dir, input_file=input_f)
         self.queue.task_done()

if __name__ == "__main__":
   # Build argument parser.
   parser = argparse.ArgumentParser(
         prog = 'filterGemSam.py',
         description = """ filter the gemsam file
         - sam files obtained by gem-2-sam
         - n: number of cores for multiprocessing
         !NOTE this program requires a big chunk of memory, the binning dumped directly into memory to speed up the filtering!""",
         formatter_class=RawTextHelpFormatter
   )
   parser.add_argument(
         '-n',
         '--n-threads',
         metavar = 't',
         type = int,
         nargs = '?',
         default = 1,
         help = 'number of threads to run in paralel'
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
 
   # filter with threads
   lock = multiprocessing.Lock()
   queue = multiprocessing.JoinableQueue(-1) #no.limit for queue
   manager = multiprocessing.Manager()


   # Fill in task queue while threads are running.
   for f_task in glob.glob(args.input_file):
      queue.put_nowait(f_task)

   # Create and start filtering in threads.
   threads = [filterThread(lock, queue, args.output_dir) for i in range(args.n_threads)]
   for thread in threads:
      thread.start()
   queue.join()

