#!/usr/bin/env python
            
##usage: ./long_branch_finder2.py <tabular treelist file> <number of stdevs/mads as cut-off>  > outfile
## JPM edited version. Calculates median absolute deviations and/or std dev of tree's branch lengths *without* QUERY sequences.
#import modules
import sys, os, numpy, re
from statsmodels import robust #JPM EDITS

def read(filename):
  f = open(filename)
  lines = f.readlines()
  for eachline in lines:
    line = eachline.split('\t')
    gene = line[0]
    d1 = {}
    d1[gene] = line[1]   #matches genename with its tree
    treetips = re.findall('[a-zA-Z0-9]+(?:_[a-zA-Z0-9]+)?\w*:\d+\.\d+', line[1]) #should be more flexible in recognizing speciesnames in trees
    # treetips = re.findall('[a-zA-Z0-9]+(?:_[a-zA-Z0-9]+)?:\d+\.\d+', line[1])  # makes a list of items like 'spname:bl'
    #treetips = re.findall('[A-Z][a-z]+_[a-z]+:\d+\.\d+', line[1])  # makes a list of items like 'spname:bl'

    pat = re.compile('[a-zA-Z0-9]+(?:_[a-zA-Z0-9]+)?\w*:\d+\.\d+') #JPM EDITS
    treetips2 = filter(lambda x: not x.startswith('QUERY'), pat.findall(line[1])) #JPM EDITS

    d2 = {}
    for i in treetips:
      spbl = i.split(':')
      d2.update({spbl[1] : spbl[0]})  #creates link betwn taxon and its BL
    tipbl = re.findall('\d+\.\d+', str(treetips))

    tipbl2 = re.findall('\d+\.\d+', str(treetips2))
    mad = robust.mad([float(i) for i in tipbl2], c=1) #JPM: Calculates median absolute deviation of non-normalized branch lengths without QUERY sequences
    std = numpy.std([float(i) for i in tipbl2]) #JPM: Calculates std dev of mean branch lengths without QUERY sequences
    #std = numpy.std([float(i) for i in tipbl]) #Original: Calculates std dev of mean branch lengths with QUERY sequences

    #print "Median Absolute Deviation of BL for", gene + ":", mad        #just
    #print                                                               #testing
    #print "Standard Deviation of the Mean of BL for", gene + ":", std   #some
    #print                                                               #stuff

    numstd = int(sys.argv[2])*mad #JPM: Change mad to std to switch back to standard deviations
    for i in tipbl:
      if float(i) > float(numstd):
#        print d2[str(i)] + '\t' + gene  + '\t' + i  
        print d2[str(i)] + '\t' + gene    

  f.close()


def main():
  read(sys.argv[1])

if __name__ == '__main__':
  main()


