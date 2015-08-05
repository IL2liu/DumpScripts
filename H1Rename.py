#!/usr/bin/python
#date of birth : 140912

import os

def H1rename(annot_file):
   """
   Rename File with adding the reference assembly 
   in front separated by '-'
   """

   annot_file = open(annot_file)
   dict_file = {}
   ref_name = ''
   file_name = ''
   for row in annot_file:
      if row[0] == '#':
         continue
      row = row.rstrip()
      row = row.split(',')
      ref_name = row[5]
      file_name = row[10]
      dict_file[file_name] = ref_name
   
   # add extra file name for the failed downloaded file H1-872b.bam
   dict_file['H1-872b.bam'] = "None"

   for filename in os.listdir("Input_Data/."):
      if filename in dict_file.keys():
         new_name = dict_file[filename] + '-' + filename
         os.rename(os.path.join("Input_Data",filename), os.path.join("Input_Data",new_name))
      continue


if __name__ == "__main__":
   H1rename("H1-hESC_annotation.txt")

      
