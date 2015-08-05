#!/usr/bin/python                        
#date of birth: 050912
#date of modification: 151012
 # modified for fasq file

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

   """ #for more than 1 cell
   cells = []
   for cell in cells_selection:
      cell = cell.string.strip()[:2] #get the first two letters of the cell name
      cells.append(cell)
   cells = '|'.join(cells)
   """
   
   search_title = "#T" + str(soup.find("title").string.strip()) +\
            "|Searched cells:" + cell + "\n"

   #header name for each column's record
   header_records = "#pi,lab,exp,view,replicate_no,assembly_map,ucsc_no,geo_no,add_info,link,file_name\n"


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
      ass_map = col[6].string.strip()
      ucsc_no = col[8].string.strip()
      geo_no = col[11].string.strip()
      add_info = col[15].string.strip()[:-1].replace(';','|') # remove the last ';'
      
      #get the download links
      re_link = re.compile('window.location=\'([a-zA-Z0-9\.-_]*)')
      link = re_link.findall(str(col[0]))

      #to make link pretty, as it should for gentlemen
      link = str(link).strip('[]').strip().replace("'", "")

      #provide the file_name
      #use the first two letters of the cell's name
      file_type = col[14].string.strip()
      file_name = '%s-%s-%03d.%s' %(ass_map, cell[:2], counter_row, file_type)
      counter_row += 1       

      record = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' %\
               (pi, lab, exp, view, replicate, ass_map, ucsc_no,
                geo_no, add_info, link, file_name)
      record = str(record).replace("&nbsp;", "None")
      
      #select only the alignments files
      """
      if view == "Alignments" and exp == "ChIP-seq":
         records.append(record)
      else:
         continue
      """
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
                  
def RetrieveLink(annot_file, suff_dir = None):

   annot_file = open(annot_file, "r")
   
   for item in annot_file:
      item = item.rstrip()
      #skipping the comments
      if item[:2] == '#T':
         #get the full searched cells' name for dir's name
         cell_name = re.search('cells:.+', item).group(0)
         cell_name = suff_dir + cell_name
         #create the dir with cell's name
         if not os.path.exists(cell_name):
            os.makedirs(cell_name)
         continue
      elif item[0] == '#':
         continue
      else:
         download_url = item.split(",")[-2]
         file_name = item.split(",")[-1]
         
         open_url = urllib2.urlopen(download_url)
         f_os = str(os.path.join(cell_name, file_name))
         with open(f_os, 'w') as f:
            #progress bar, print downloading and status info
            meta = open_url.info()
            file_size = int(meta.getheaders("Content-Length")[0])
            print "Currently downloading: %s Bytes: %s" %\
                  (file_name, file_size)

            file_size_dl = 0
            block_sz = 8192
            while True:
               buffer = open_url.read(block_sz)
               if not buffer:
                  break
               file_size_dl += len(buffer)
               f.write(buffer)
               status = r"%10d  [%3.2f%%]" %\
                        (file_size_dl, file_size_dl * 100. / file_size)
               status = status + chr(8)*(len(status)+1)
               print status
               #sys.stderr.write(status + "\n")

if __name__ == '__main__':

   #getting human H1-hESC fasq_url
   fasq_url = "http://encodeproject.org/cgi-bin/hgFileSearch?hgsid=304019069&db=hg19&hgt_tsDelRow=&hgt_tsAddRow=&tsName=&tsDescr=&tsGroup=Any&fsFileType=fastq&hgt_mdbVar1=cell&hgt_mdbVal1=H1-hESC&hgt_mdbVar2=antibody&hgt_mdbVal2=Any&hgfs_Search=search"
   
   #parse the searchtable
   fasq_tableRecords = parseSearchTable(fasq_url)
   
   #write the annotation file for H1-hESC
   writeAnnot(fasq_tableRecords, 'Fasq-H1-hESC')

   #get the files
   #run the script with 2> log.txt
   fasq_download_file = RetrieveLink("Fasq-H1-hESC_annotation.txt", "fasq-")
