import os
import optparse
import subprocess
from multiprocessing import Pool

directory = ""
results = "results.data"
extension = ""
aligned_extension = ".tab"
datatype = ""

perlpath = "/home/galaxy-dist/tools/osiris/tree-manipulation/"

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

#def toData(text, name):
#    name = name.replace("fasta", "") #file name has fasta when fasta file called
#    text = name.replace(".fs.tre", "") + "\t" + text.replace(" " , "")
#    return text


def toData(text, name):
    text = text.split('\n') 
    result = ''
    for line in text:
        if '\t' in line:
        	line = line.replace("./data/","") + "\n"
        result += line
    return result  # Index past the first newline char

def LB_pruner(input):
    file_name = directory + os.sep + input
    popen = subprocess.Popen(['perl', perlpath+'LB_prunerG.pl', file_name, indata, file_name + aligned_extension])
    popen.wait()

class Sequence:
    def __init__(self, string):
        lis = string.split()
        self.name = lis[0]
        self.tree = lis[1]
        self.string = string

    def printFASTA(self):
        return self.tree + '\n'

def saveMulti(tabFile):
    with open(tabFile) as f:
        for line in f:
            seq = Sequence(line)
            with open(directory + os.sep + seq.name + extension, "a") as p:
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
        '-m', '--mult',
        dest='datatype',
        action='store',
        type='string',
        help='Multiplier')

    options, args = parser.parse_args()

    global directory
    global indata
    inputFile = unescape(options.input)
    directory = unescape(options.path) + os.sep + "data"
    indata = unescape(options.datatype)

    os.mkdir(directory)

    if isTabular(inputFile):
        saveMulti(inputFile)
    else:
        saveSingle(inputFile)

    pool = Pool()
    list_of_files = [file for file in os.listdir(directory) if file.lower().endswith(extension)]
    pool.map(LB_pruner, list_of_files)

    result = [file for file in os.listdir(directory) if file.lower().endswith(aligned_extension)]
    with open(directory + os.sep + results, "a") as f:
        for file in result:
            with open(directory + os.sep + file, "r") as r:
                f.write(toData(r.read(),file))

if __name__ == '__main__':
    main()

