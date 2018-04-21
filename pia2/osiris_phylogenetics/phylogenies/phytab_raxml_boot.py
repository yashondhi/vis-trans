#!/usr/bin/env python

##  phytab_raxml_boot
##
##  This tool runs RAxML on each partition/gene in a PHYTAB file.
##  If N = # of nodes requested in job runner, then N RAxML jobs will run simultaneously.  Make sure that the
##  number of processors ('ppn') in the job runner matches the 'numthreads' value set in command line arguement.
## 
##  Usage:   ./phytab_raxml.py -i <phytabinput> -e <model> -f <modelfile> -T 4
##  example: ./phytab_raxml.py -i myphytab.txt -e PROTGAMMAWAG -f None -T 4
##  or:      ./phytab_raxml.py -i myphtab.txt -e None -f modelsforeachpartition.txt -T 4

import optparse
import os
import subprocess
import multiprocessing

RESULTS_DIR = 'results'
RESULTS_FILE = 'results.phy'
RAXML_PREFIX = 'RAxML_bipartitionsBranchLabels.'

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


class Species:
    def __init__(self, string):
        lis = string.split('\t')
        # print lis
        self.species = lis[0]
        self.gene = lis[1]
        self.name = lis[2]
        self.sequence = lis[3]

    def toString(self):
        return self.species + '\t' + self.sequence


class Gene:
    def __init__(self, name):
        self.name = name
        self.count = 0
        self.length = 0
        self.species = []

    def output(self):
        file_name = self.name + ".phy"
        location = RESULTS_DIR + os.sep + file_name
        with open(location, 'w') as f:
            f.write(str(self.count) + '\t' + str(self.length) + '\n')
            for s in self.species:
                f.write(s.toString())
        return file_name

    def add(self, species):
        if species.name == "":
            return
        self.species.append(species)
        self.count += 1
        if self.length == 0:
            self.length = len(species.sequence) - 1


def output_species(species):
    file_name = species.gene + ".phy"
    location = RESULTS_DIR + os.sep + file_name
    with open(location, 'a') as f:
        f.write(species.toString())
    return file_name


def process_phytab(input):
    files = set()
    genes = dict()
    with open(input) as f:
        for line in f:
            if len(line) < 4:
                continue
            species = Species(line)
            if species.gene in genes:
                genes[species.gene].add(species)
            else:
                gene = Gene(species.gene)
                gene.add(species)
                genes[gene.name] = gene
    for k, gene in genes.iteritems():
        files.add(gene.output())
    return files


def runRaxml(list_of_files, evo, evoDict, NUMTHREADS):
  count = 0
  list_of_files = sorted(list_of_files)
  for gene_file in list_of_files:
    count+=1            
    if gene_file.split(".")[0] in evoDict:
      newEvo = evoDict[gene_file.split(".")[0]]
    else:
      newEvo = evo
    file_name = RESULTS_DIR + os.sep + gene_file
## RAxML notes:  ##  
# to run parsimony trees:
#          popen = subprocess.Popen(['raxmlHPC-PTHREADS', '-T', cpu_count,'-f', 'd', '-s', file_name,'-y', '-m',  newEvo, '-n', gene_file[:-4]+'.tre', '-p', '34'])
# to run likelihood trees:
#          popen = subprocess.Popen(['raxmlHPC-PTHREADS', "-T", cpu_count, "-s", file_name, '-m', newEvo, '-n', gene_file[:-4], '-p', '34'])
# to run likelihood trees using starting tree:
#    popen = subprocess.Popen(['raxmlHPC-PTHREADS', '-T', cpu_count, '-f', 'e','-s', file_name, '-m', newEvo, '-n', gene_file[:-4], '-t', ptre])   

## run trees as simultaneously as possible--but wait until last tree completes to resume script
    raxml_cmd = ['raxmlHPC-PTHREADS-SSE3', '-T', NUMTHREADS, '-f', 'a',
                 '-s', file_name, '-m', newEvo, '-n', gene_file[:-4],
                 '-p', '1234567', '-x', '1234567', '-N', '100', '-K', 'MK']
    if count == len(list_of_files):
      run = subprocess.Popen(raxml_cmd)
      run.wait()
    else:
      run = subprocess.Popen(raxml_cmd)
      run.communicate()[0]


def readEfile(efile):
    evoDict = {}
    with open(efile, "r") as f:
        for line in f:
            pair = line.split("\t")
            evoDict[pair[0].strip()] = pair[1].strip()
    return evoDict

def main():
    usage = """%prog [options]
options (listed below) default to 'None' if omitted
    """
    parser = optparse.OptionParser(usage=usage)

    parser.add_option(
        '-i', '--in',
        dest='input',
        action='store',
        type='string',
        metavar="FILE",
        help='Name of input data.')

    parser.add_option(
        '-e', '--evo',
        dest='evo',
        action='store',
        type='string',
        metavar="EVO",
        help='Evolution model.')

    parser.add_option(
        '-f', '--evo-file',
        dest='efile',
        action='store',
        type='string',
        metavar="EVO_FILE",
        help='Evolution model file. Format is gene_name [tab] evolution_model.')

    parser.add_option(
        '-T', '--numthreads',
        dest='numthreads',
        action='store',
        type='int',
        metavar="NUMT",
        help='Specify number of threads.')
        
    options, args = parser.parse_args()

    os.mkdir(RESULTS_DIR)

    list_of_species_files = process_phytab(unescape(options.input))

    try:
        evoDict = readEfile(unescape(options.efile))
    except IOError:
        print "No sequence model file provide...using", unescape(options.evo), "as the model" 
        evoDict = {}
   
    runRaxml(list_of_species_files, unescape(options.evo), evoDict, str(options.numthreads))
    
    result = [file for file in os.listdir('./') if file.startswith(RAXML_PREFIX)]
    result = sorted(result)
    with open(RESULTS_DIR + os.sep + RESULTS_FILE, "w") as f:
        for file in result:
            with open(file, "r") as r:
                f.write(file[len(RAXML_PREFIX):] + '\t' + r.read())
        f.close()

if __name__ == '__main__':
    main()

