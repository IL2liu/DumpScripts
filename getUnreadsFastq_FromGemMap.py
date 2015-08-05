#!/usr/bin/env python
# -*- coding: utf-8 -*-

#date of created: 261212

############################################
# get the unmapped reads from the gem-map  #
############################################

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

from gzopen import gzopen

def get_Unmap(input_file, output_dir='.'):
   """ get the unmapped reads from the gemMap output """

   try:
      with gzopen(input_file) as in_f:
         gc.disable()
         reads_unmapped = []
         for line in in_f:
            items= line.rstrip().split('\t')
            if items[3] =="0:0":
               # collect all the fasqt elements according to the fastq format (ENCODE)
               # 1. the ids of the reads
               # 2. the read seq
               # 3. the strand + the ids
               # 4. the quality code
               # the unmapped regions of the gemMap are given '+' strand(assumption)
               fastq_items = []
               fastq_items.extend(["@"+items[0], items[1], "+"+items[0], items[2]])
               reads_unmapped.append(fastq_items)
         gc.enable()
   except:
      sys.stderr.write('file error:%s'%(os.path.splitext(input_file)[0]))
      raise

   # dump for the output...
   if not os.path.exists(output_dir):
      try:
         os.makedirs(output_dir)
      except OSError as exception:
         if exception.errno != errno.EEXIST:
            raise
   head, tail = os.path.split(input_file)
   base = os.path.splitext(tail)[0]
   output_fname_unmap ="unmapped-%s.fastq" %(base)
   output_file_unmap = str(os.path.join(output_dir, output_fname_unmap))

   with open(output_file_unmap, 'w') as output_f_unmap:
      for line in reads_unmapped:
         for items in line:
            output_f_unmap.write(items+'\n')

class getUnmapThread (multiprocessing.Process):
   """ parse the unmapped reads in threads """
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

         get_Unmap(output_dir=self.output_dir, input_file=input_f)
         self.queue.task_done()


if __name__ == "__main__":
   # Build argument parser.
   parser = argparse.ArgumentParser(
         prog = 'getUnreadsFastq_FromGemMap.py',
         description = """get unmapped reads 
         !NOTE this program requires a big chunk of memory, the parsed unmapped reads dumped directly into memory to speed up the process!""",
         formatter_class=RawTextHelpFormatter
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

   # parse the unmapped reads with threads 
   lock = multiprocessing.Lock()
   queue = multiprocessing.JoinableQueue(-1) #no.limit for queue 
   manager = multiprocessing.Manager()
 
   # Fill in task queue while threads are running. 
   for f_task in glob.glob(args.input_file):
      queue.put_nowait(f_task)
 
   # Create and start binning in threads. 
   threads = [getUnmapThread(lock, queue, args.output_dir) for i in range(args.n_threads)]
   for thread in threads:
      thread.start()
   queue.join()
