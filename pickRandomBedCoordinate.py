#!/usr/bin/env python
# -*- coding:utf-8 -*-

#date of birth:200912

#######################################
# Randomly picked the bed coordinates #
# default to picked 1000, line_nos    #
#######################################

import sys
import random
import argparse
import os

def pickBedCoordinateRandomly(bed_file, line_nos=1000, dir_store="Output_Data"):
   """
   generate bed file with randomly chosen coordinates(randomly picked the line)
   and return line_list
   This returned line_list is used to get the bed coordinates
   from another file for sequence comparison
   !!! dir_store must be similar with the bash script running this py scripy
   """

   #replace the .bed with -random.bed as the output file names
   randBed_file = "randomlyPicked-" + os.path.basename(bed_file)

   #read all the lines
   with open(bed_file) as f:
      all_lines = f.readlines()

   if type(line_nos) is int:
      line_nos = random.sample(xrange(len(all_lines)), line_nos)

   sample = [all_lines[i] for i in line_nos]

   with open(os.path.join(dir_store,randBed_file) , "w") as a:
      for line in sample:
         a.write(line)

   return line_nos


if __name__ == "__main__":


   if (len(sys.argv[0:]) < 3) or sys.argv[1]=="-h":
      print """ usage: pickRandomBedCoordinate.py [f_after] [f_before ]
        f_after            please add your bedfile after lifting over
        f_before           please add bedfile before lifting over """
      exit(-1)

   #pick randomly the bed coordinate of the file after lifted over
   line_list = pickBedCoordinateRandomly(sys.argv[1])
   
   #pick the bed coordinate of the file before lifting over
   pickBedCoordinateRandomly(sys.argv[2], line_list)
