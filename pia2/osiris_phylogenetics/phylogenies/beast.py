#!/usr/bin/env python
"""
This program makes a new copy of a BEAST XML config file and changes the name of the log and tree files.
The names have to be changed since Galaxy does not support dynamic output file names.
The new XML file name must not be the same as the original!!!
Script prints unique filename.


Usage: python beast.py XMLFILE SAVE_DIRECTORY
"""

import sys
import time
from xml.etree.ElementTree import ElementTree

beastDOM = ElementTree()
beastDOM.parse(sys.argv[1])
logs = beastDOM.findall('mcmc/log')
for log in logs:
    if log.get('id', False) == "fileLog":
        log.set('fileName', 'data.log')

logs = beastDOM.findall('mcmc/logTree')
for log in logs:
    if log.get('id', False) == "treeFileLog":
        log.set('fileName', 'data.trees')

if len(sys.argv) > 2:
    directory = sys.argv[2]
else:
    directory = ""

filename = directory + "/" + str(time.time()) + '.xml'
beastDOM.write(filename)
print filename
