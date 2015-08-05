#!/usr/bin/python                        
# -*- coding:utf-8 -*-

#######################################
# Parse ENCODE download search table  #
#######################################

import sys
import urllib2
from BeautifulSoup import BeautifulSoup
import re
import csv
import os
      
def parseSearchTable(url):
   open_html = urllib2.urlopen(url)
   html = open_html.read()
   soup = BeautifulSoup(html)

   #get the selected cells during the search
   selection = soup.find("select", {"name":"hgt_mdbVal1"})
   cells_selection = selection.findAll("option",
                     {"selected":"SELECTED"})
   #get the cell
   cell = cells_selection[0].string.strip()

   search_title = "#T" + str(soup.find("title").string.strip()) +\
            "|Searched cells:" + cell + "\n"

   #header name for each column's record
   header_records = "#pi,lab,exp,view,replicate_no,ucsc_no,geo_no,add_info,link,file_name\n"

   #get the selected tables
   tableBody = soup.find("tbody", {"class":"sortable sorting"})

   #collect data start with the search title
   records = [search_title, header_records]

   counter_row = 1
   for row in tableBody.findAll("tr"):
      col = row.findAll("td")
      
      #fields stored
      pi = col[1].string.strip()
      lab = col[2].string.strip()
      exp = col[3].string.strip()
      view = col[4].string.strip()
      replicate = col[5].string.strip()
      ucsc_no = col[7].string.strip()
      geo_no = col[10].string.strip()
      add_info = col[13].string.strip()[:-1].replace(';','|') # remove the last ';'
      add_info = add_info.replace(","," ")
      
      #get the download links
      re_link = re.compile('window.location=\'([a-zA-Z0-9\.-_]*)')
      link = re_link.findall(str(col[0]))

      #to make link pretty, as it should for gentlemen
      link = str(link).strip('[]').strip().replace("'", "")

      #provide the file_name
      #use the first two letters of the cell's name
      file_type = col[12].string.strip()
      file_name = '%s-%03d.%s' %(cell[:2], counter_row, file_type)
      counter_row += 1       

      record = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' %\
               (pi, lab, exp, view, replicate, ucsc_no,
                geo_no, add_info, link, file_name)
      record = str(record).replace("&nbsp;", "None")
      
      #select only the alignments files
      if view == "Raw data" and exp == "ChIP-seq":
         records.append(record)
      else:
         continue
   return records

def writeAnnot(records, annote_name):

   #write output file as text and csv files
   with open(annote_name + '_annotation.txt', "w") as f:
      for item in records:
         f.write(item)
                 
   f_csv =  csv.writer(open(annote_name + "_annotation.csv", "wb"))
   f_txt = csv.reader(open(annote_name + '_annotation.txt', "rb"), delimiter = ",")
   f_csv.writerows(f_txt)
                  
if __name__ == '__main__':

   #parse the searchtable
   fasq_tableRecords = parseSearchTable(sys.argv[1])
   
   #write the annotation file for H1-hESC
   writeAnnot(fasq_tableRecords, sys.argv[2])
