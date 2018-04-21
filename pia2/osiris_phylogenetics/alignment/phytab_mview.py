#!/usr/bin/env python
## usage: ./phytab_mview.py -i <phytabinput> -d <protein|dna> 
## splits up an aligned phytab file containing multiple genes into
## individual files to run mview

import sys, os, os.path, tempfile, shutil, re, shlex, subprocess
import optparse
from multiprocessing import Pool

#define some variables to call later:

directory = ""
extension = ".fs"
html_header = """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
<HTML>
<HEAD>
<TITLE></TITLE>
</HEAD>
<BODY BGCOLOR='white' TEXT='black' LINK='blue' ALINK='red' VLINK='purple'>
<H1>PHYTAB MVIEW ALIGNMENT VIEWER</H1>
<PRE>Select from below to view aligned sequence as HTML (left) or FASTA (right) in browser.
</PRE>
<table border="1" bordercolor="#000000" style="background-color:#FFFFFF" width="300" cellpadding="3" cellspacing="0">
	<tr>
		<td>mview HTML</td>
		<!--<td>FASTA</td>-->
	</tr>"""
html_close =  """
<P><SMALL><A HREF="http://bio-mview.sourceforge.net">MView</A> </SMALL><BR>
</BODY>
</HTML>"""	

#define some functions to call in 'main':
#    first, sanitize problematic characters
def unescape(string):
  mapped_chars = {
        '>': '__gt__',
        '<': '__lt__',
        "'": '__sq__',
        '"': '__dq__',
        '[': '__ob__',
        ']': '__cb__',
        '{': '__oc__',
        '}': '__cc__',
        '@': '__at__',
        '\n': '__cn__',
        '\r': '__cr__',
        '\t': '__tc__',
        '#': '__pd__'
        }

  for key, value in mapped_chars.iteritems():
    string = string.replace(value, key)

  return string
#  next, define tabular --> fasta conversion
class Sequence:            
  def __init__(self, string):
    lis = string.split()
    self.species = lis[0]
    self.family = lis[1]
    self.name = lis[2]
    self.header = ' '.join(lis[:-1])
    self.sequence = lis[-1]
    self.string = string

  def printFASTA(self):
    return '> ' + self.header + '\n' + self.sequence + '\n'

#  then define function to apply preceding conversion method to all genes
#  (creates separate file for each gene)
def saveMulti(tabFile):
  with open(tabFile) as f:
    for line in f:
      seq = Sequence(line)
      with open(seq.family + extension, "a") as p:
        p.write(seq.printFASTA())
                
#subroutine to write main HTML output containing valid urls to mview htmls
def resultsto_output_html(html_mainoutput,basepath):
  htmllist = [f for f in os.listdir(basepath) if 'html' in f]
  sortedhtmllist = sorted(htmllist)
  html = open(html_mainoutput, 'w')
  html.write(html_header)
  for f in sortedhtmllist:
    f_path = os.path.join(basepath,f)
    htmllink = '<tr><td><a href="' + f + '">' + f + '</a></td>\n' 
    html.write(htmllink)
  html.write(html_close)
  html.close()

def main():
#the command line arguments from the xml:
  """
           ##params for galaxy wrapper
           $input 
           $dna  
           $output 
           "$output.extra_files_path"  #save the htmlfiles here
  """ 
  inputphytabfile = sys.argv[1]
  dnaorprotein = sys.argv[2]
  output = sys.argv[3]
  extra_files_path = sys.argv[4]
  
  inputFile = unescape(inputphytabfile)
  ##make the fasta files
  saveMulti(inputFile) 

  #prepare to put mview htmls into valid path

  if not os.path.isdir(extra_files_path):  #make filepath for alns to go with galaxy info
      os.makedirs(extra_files_path)    
  
  # execute mview on each fasta, storing in extra_files_path as <gene_aln>.html
  list_of_fastafiles = [f for f in os.listdir(os.getcwd()) if 'fs' in f]
  sortedfileorder = sorted(list_of_fastafiles)
  for gene_aln in sortedfileorder:
    result_htmlfile = gene_aln + '.html'
    result_path = os.path.join(extra_files_path,result_htmlfile) #puts the htmls in permanent Galaxy directory
    if dnaorprotein is 'dna':
      cmd = subprocess.Popen(['mview','-in','pearson','-DNA','-bold','-coloring','group','-html','head', gene_aln],stdout=subprocess.PIPE)
    else:
      cmd = subprocess.Popen(['mview','-in','pearson','-bold','-coloring','group','-html','head', gene_aln],stdout=subprocess.PIPE)
    cmd.wait()  
    out =  cmd.communicate()[0]
     
    with open(result_path, 'wb') as fileout:
      fileout.write(out)
    ##now have # of gene htmls in extra_files_path/
    
  #write main html output  
  resultsto_output_html(output,extra_files_path)


if __name__ == '__main__':
    main()

