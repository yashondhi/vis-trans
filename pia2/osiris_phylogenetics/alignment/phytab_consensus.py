#!/usr/bin/env python

import os
import optparse
import subprocess

import decimal

from Bio import SeqIO
from Bio.Align import AlignInfo
from Bio import AlignIO


directory = ""
extension = ""


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

def consensus(input,percent,ambiguity):
    file_name = directory + os.sep + input
#    if os.path.getsize(file_name) > 0:	#check that there is sequence in the file
    alignment = AlignIO.read(file_name, "fasta")
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.dumb_consensus(threshold=percent, ambiguous=ambiguity, require_multiple=1)
    conseq = consensus.rstrip('\n')
    return(conseq)
 
class Sequence:
    def __init__(self, string):
        lis = string.split()
        self.species = lis[0]
        self.family = lis[1]
        self.name = lis[2]
        self.header = ' '.join(lis[:-1])
        try:
            self.sequence = lis[3]
        except IndexError:
            self.sequence = 'empty'
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
        '-a', '--ambiguity',
        dest='ambiguity',
        action='store',
        type='string',
        help='Character for ambiguity.')

    parser.add_option(
        '-p', '--percent',
        dest='percent',
        action='store',
        type='int',
        help='Minimum Percent to contribute to consensus.')

    options, args = parser.parse_args()

    global directory
    inputFile = unescape(options.input)
    directory = unescape(options.path) + os.sep + "data"
#    decimalpercent = float(options.percent)
    decimalpercent = float(0.1)
    ambiguity = options.ambiguity

    os.mkdir(directory)

    if isTabular(inputFile):
        saveMulti(inputFile)
    else:
        saveSingle(inputFile)

    list_of_files = [file for file in os.listdir(directory) if file.lower().endswith(extension)]
    for gene in list_of_files : 
        current_consensus = consensus(gene, decimalpercent, ambiguity)
        print(gene + "\t" + current_consensus + "\t" + str(decimalpercent) + " " + str(options.percent))

if __name__ == '__main__':
    main()
