#!/usr/bin/env python

import os
import optparse
import subprocess
from multiprocessing import Pool
from functools import partial

directory = ""
results = "results.data"
extension = ".fs"
slim_extension = ".slim"
jarpath='/home/galaxy/pkgs/SSR/SimilarSequenceRemover.jar'

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


def isTabular(file):
    with open(file) as f:
        for line in f:
            if line[0] == '>':
                return False
    return True


def toData(text):
	text = text.split('\n')
	result = ''
	for line in text:
	    if '>' in line:
	        line = '\n' + line.replace('> ', "") + '\t'
	    line = line.replace(" ", "\t")
	    result += line
	return result[1:]  # Index past the first newline char

def toDataSingle(text):
	text = text.split('\n')
	result = ''
	for line in text:
		line = line + '\n'
 		result += line
	return result[1:]  # Index past the first newline char

def runSimilarSequenceRemover(input):
    input_name = directory + os.sep + input
    outfasta = input_name + slim_extension
    call = subprocess.call(['java','-jar', jarpath, '-i', input_name, '-o', outfasta, 's=' + percentage, '-h', alignallornot ])

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


def saveMulti(tabFile):
    with open(tabFile) as f:
        for line in f:
            seq = Sequence(line)
            with open(directory + os.sep + seq.family + extension, "a") as p:
                p.write(seq.printFASTA())


def saveSingle(fastaFile):
    with open(fastaFile) as f:
        for line in f:
            with open(directory + os.sep + "fasta" + extension, "a") as p:
                p.write(line)


def main():
    usage = """%prog [options]
options (listed below) default to 'None' if omitted
    """
    parser = optparse.OptionParser(usage=usage)

    parser.add_option(
        '-d', '--directory',
        metavar="PATH",
        dest='path',
        default='.',
        help='Path to working directory.')

    parser.add_option(
        '-i', '--in',
        dest='input',
        action='store',
        type='string',
        metavar="FILE",
        help='Name of input data.')
    
    parser.add_option(
        '-a', '--alignall',
        dest='alignall',
        action='store',
        type='string', help='t or f. t  only aligns first 100 sites instead of entire seq')
        
    parser.add_option(
        '-s', '--similarity',
        dest='similarity',
        action='store',
        type='float',
        help='Percentage similarity as cutoff (eg 0.99).')

    options, args = parser.parse_args()

    global directory
    global percentage
    global alignallornot
    inputFile = unescape(options.input)
    directory = unescape(options.path) + os.sep + "data"
    percentage = str(options.similarity)
    alignallornot = options.alignall
    
    os.mkdir(directory)

    if isTabular(inputFile):
        saveMulti(inputFile)
    else:
        saveSingle(inputFile)
    
    pool = Pool()
    list_of_files = [file for file in os.listdir(directory) if file.lower().endswith(extension)]
    list_of_files = sorted(list_of_files)
    
    pool.map(runSimilarSequenceRemover, list_of_files)
   
    result = [file for file in os.listdir(directory) if file.lower().endswith(slim_extension)]
    if isTabular(inputFile):
	with open(directory + os.sep + results, "a") as f:
	    for file in result:
	        with open(directory + os.sep + file, "r") as r:
	            f.write(toData(r.read()) + "\n")
    else:
	with open(directory + os.sep + results, "a") as f:
	    for file in result:
	        with open(directory + os.sep + file, "r") as r:
	            f.write(r.read() +'\n')

if __name__ == '__main__':
    main()


