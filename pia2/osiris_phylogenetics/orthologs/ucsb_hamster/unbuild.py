#!/usr/bin/env python
'''
Usage: unbuild.py outputDirectoryName hmmFileWithTabs
This script takes a tab sperated, newline escaped HMM models file then splits the models up and makes a file list.
The output files are put into the directory outputDirectoryName.
'''
from sys import argv, exit
import os


def parse(line):
    line = line.partition('\t')
    name = line[0] + ".hmm"
    data = line[2].replace("\\n", "\n")
    return name, data


def writeFile(location, data):
    try:
        with open(location, 'w') as f:
            f.write(data)
    except:
        pass


def createDir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)


if __name__ == "__main__":
    if len(argv) != 3:
        print "Usage: unbuild.py outputDirectoryName hmmFileWithTabs"
        exit()
    dir = argv[1]
    createDir(dir)
    
    name_list = [] #  Stores new file names

    #  Read data and seperate into files.
    with open(argv[2], 'r') as f:
        for line in f:
            name, data = parse(line)
            if data.strip() == "":  # no data? skip it.
                continue
            name_list.append(name)
            name = dir + os.sep + name
            writeFile(name, data)

    #  Write new file names to a file.
    name_list_location = dir + os.sep + "hmmlist.txt"
    with open(name_list_location, 'w') as f:
        for name in name_list:
            f.write(name + "\n")
