#!/usr/bin/python
# -*- coding:utf-8 -*-

import os
import glob
from collections import defaultdict

def binit(dirBam, window_size=None):
   """
   binit the alignment files(.bam)
   all the bam files bin-indexed starts with 1-max'bin index
   bin's index = int(alignment position/window_size)
   max's bin index = max(alignment position/window_size)

   """
   
   #get the bam files
   path = '../%s'%(dirBam) #dirBam: where the bam files are located

   #get all the .bam files
   for bamfile in glob.glob(os.path.join(path, '*.bam')):
      
      #name of the bam files without dir names
      bam_txt = bamfile.split('/')[-1]

      try:
      #create txt file from bam file using samtools 
         os.system('samtools view %s > %s.txt'%(bamfile,bam_txt))
      except:
         print "Unexpected error:", sys.exc_info()[0]
         raise

      #binning starts ....
      f_binName = bam_txt +'_bin'
      with open(f_binName,'w') as f_bin:
         #store chromosome, and bin's index as keys
         #read coverage as values
         freq_table = defaultdict(int)
         
         f = open(bam_txt+'.txt','r')
         for line in f:
            items = line.split()
            freq_table[(items[2],int(items[3])/window_size)] += 1

         #collect the chromosomes without duplicates
         chroms = set([chrom for (chrom,bin) in freq_table.keys()])
         for chrom in chroms:
            maxbin = max([bin for (chrom__,bin) in freq_table.keys() if chrom__ == chrom])

            for bin in range(maxbin):
               f_bin.write("%s\t%s\t%s\t%d\n" % \
                     (chrom, 1+bin*window_size, (1+bin)*window_size,
                        freq_table[(chrom,bin)]))
         #zip the file
         #os.system('gzip %s' %(f_binName))

if __name__ == "__main__":
   #self-test code
   binit("cells", 300)


