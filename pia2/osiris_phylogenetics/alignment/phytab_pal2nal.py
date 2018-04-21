#!/usr/bin/env python

"""
pal2nal.py
----------
Runs pal2nal.pl on GENENAME protein/DNA file pairings. Each file with protein
    information must have a corresponding file with DNA information for its
    gene. Depends on the Splitter class from phytab_splitter library which may 
    be installed form PyPI. Designed for use with the Galaxy
    bioinformatics platform.
Place pal2nal.pl in the same directory as the tool. Edit the CDIR line so that
    it points to the directory where the tool is installed.
"""

import os
import sys
import zipfile
import shutil
from splitter import Splitter

__author__ = 'William Chen'
__license__ = 'BSD 2-clause'
__status__ = 'In development'

# declare working directory
CDIR = '/home/galaxy2/galaxy-dist/tools/osiris_phylogenetics/alignment/'

# get inputs for pal2nal.pl
protFile = sys.argv[1] # Phytab file with protein info
dnaFile = sys.argv[2] # Phytab file with DNA info 
outputType = sys.argv[3]
blockonly = sys.argv[4]
nogap = sys.argv[5]
nomismatch = sys.argv[6]
codontable = sys.argv[7]
outFormat = sys.argv[8]

# remove whatever is at the file location marked as the output, if exists 
try:
    os.remove(sys.argv[9])
except OSError:
    pass
# change working directory to where we want
os.chdir(CDIR)
# call Splitter from splitterforpal2nal.py and split the Phytab files into
#  FASTA files
sp = Splitter(dnaFile=dnaFile, protFile=protFile)
sp.generateFasta()
# create zip file to collect pal2nal.pl output
zip = zipfile.ZipFile('myzip.zip', 'w')
# run FASTA output through pal2nal.pl and write pal2nal.pl's output into a zip
#  file
# TODO: Use Popen instead of os.system
for key in sp.getProteinInfo():
    os.system('perl pal2nal.pl ' + key + 'Prot.fas ' + key + 'DNA.fas' + ' ' +
              outputType + ' ' + blockonly + ' ' + nogap + ' ' + nomismatch +
              ' ' +  codontable + ' -output ' + outFormat + ' 1> out' + key +
              ' 2>/dev/null')
    zip.write('out' + key)

# move the zip file to the output location
os.system('mv myzip.zip ' + sys.argv[9])

# print the output location for reference
print(sys.argv[9])
