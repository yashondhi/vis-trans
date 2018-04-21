#!/usr/bin/env python

"""
splitforpal2nal.py: Splits a Phytab file into FASTA files; designed for use
    with pal2nal script
"""

import optparse
import os

__author__ = 'William Chen'
__status__ = 'In development'


class Splitter:
    """
    Class structure to convert Phytab files into FASTA format
    """

    def __init__(self, dnaFile='', protFile='', minSeqSize = 0):
        """Constructor; opens input files and reads their lines"""
        self.processProteins = False
        self.processDNA = False
        self.minSequenceSize = minSeqSize
        if protFile:
            self.fileProteins = open(protFile)
            self.proteinsLines = self.fileProteins.readlines()
            self.processProteins = True;
        if dnaFile:
            self.fileDNA = open(dnaFile)
            self.DNALines = self.fileDNA.readlines()
            self.processDNA = True;

    def getSpecies(self, lines):
        """
        Gets species from each phytab line; does this by filling up a new array
          with the first element of each line
        """
        species = []
        for line in lines:
            species.append(line.split()[0])
        return species

    def getGeneNames(self, lines):
        """
        Gets gene names from each phytab line; does this by filling up a new
          array with the second element of each line
        """
        names = []
        for line in lines:
            names.append(line.split()[1])
        return names

    def getAlignedData(self, lines):
        """
        Gets the alignment data for each phytab line; does this by filling up a
          new array with the last element of each line
        """
        alignedData = []
        for line in lines:
            tokenizedLine = line.split()
            alignedData.append(tokenizedLine[len(tokenizedLine) - 1])
        return alignedData

    def getGeneId(self, lines):
        geneIds = []
        for line in lines:
            tokenizedLine = line.split()
            geneIds.append(tokenizedLine[2])
        return geneIds

    def getProteinInfo(self):
        """
        Gets the species name, gene name, and alignment data for every protein;
          for each line, puts this data into a tuple, then collects all the tuples
          in an array
        """
        allInfo = []
        info = {}
        proteinLines = self.proteinsLines
        species = self.getSpecies(proteinLines)
        names = self.getGeneNames(proteinLines)
        geneIds = self.getGeneId(proteinLines)
        data = self.getAlignedData(proteinLines)
        uniqueNames = set(names)
        for i in range(0, len(data)):
            allInfo.append((species[i], names[i], geneIds[i], data[i]))
        for uniqName in uniqueNames:
            info[uniqName] = []
            for item in allInfo:
                if item[1] == uniqName:
                    info[uniqName].append([item[0], item[2], item[3]])
        return info

    def getDNAInfo(self):
        """
        Gets the species name, gene name, and alignment data for every DNA
          sequence for each line, puts this data into a tuple, then collects all
          the tuples in an array
        """
        allInfo = []
        info = {}
        DNALines = self.DNALines
        species = self.getSpecies(DNALines)
        names = self.getGeneNames(DNALines)
        geneIds = self.getGeneId(DNALines)
        data = self.getAlignedData(DNALines)
        uniqueNames = set(names)
        for i in range(0, len(data)):
            allInfo.append((species[i], names[i], geneIds[i], data[i]))
        for uniqName in uniqueNames:
            info[uniqName] = []
            for item in allInfo:
                if item[1] == uniqName:
                    info[uniqName].append([item[0], item[2], item[3]])
        return info

    def writeFasta(self, info, dataName):
        """
        Writes out fasta files for each gene name based on info passed in; info
          can either be DNA or protein info
        """
        for key in info:
            if len(info[key]) < self.minSequenceSize:
                continue
            f = open(key + dataName + '.fas', 'w')
            for tup in info[key]:
                f.write('>' + tup[0] + key + tup[1] + '\n' + tup[2] + '\n')

    def generateFasta(self):
        """Using writeFasta(), generates the fasta files we want"""
        if self.processProteins:
            self.writeFasta(self.getProteinInfo(), 'Prot')
        if self.processDNA:
            self.writeFasta(self.getDNAInfo(), 'DNA')


def main():
    # parses input
    parser = optparse.OptionParser(description='Split phytab files for pal2nal use')
    parser.add_option('--dna', dest='dnaFileName', help='Phytab file containing DNA sequences')
    parser.add_option('--prot', dest='protFileName', help='Phytab file containing protein alignments')
    parser.add_option('--getdinfo', dest='getdinfo', help='Get DNA info')
    (options, args) = parser.parse_args()
    # initializes Splitter class
    sp = Splitter(dnaFile=options.dnaFileName, protFile=options.protFileName)
    if options.getdinfo:
        print(sp.getDNAInfo())
    # calls generateFasta() to generate the fasta files we want
    sp.generateFasta()

if __name__ == '__main__':
    main()
