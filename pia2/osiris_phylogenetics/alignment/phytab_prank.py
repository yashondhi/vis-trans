import os
import optparse
import subprocess
from multiprocessing import Pool

directory = ""
results = "results.data"
extension = ".fs"
aligned_extension = ".afa"
output_extension = ".afa.2.fas"


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
            line = '\n' + line.replace('>__XX__', "") + '\t'
        line = line.replace("__XX__", "\t")
        result += line
    return result[1:]  # Index past the first newline char


def prank(input):
    file_name = directory + os.sep + input
    popen = subprocess.Popen(['pwd'])  
    popen.wait()
    popen = subprocess.Popen(['prank', "-d=" + file_name, "-o=" + file_name + aligned_extension, "-quiet"])  
    popen.wait()

class Sequence:
    def __init__(self, string):
        lis = string.split()
        self.species = lis[0]
        self.family = lis[1]
        self.name = lis[2]
        self.header = '__XX__'.join(lis[:-1]) #prank replaces space with _ so can't join with spaces like muscle does
        self.sequence = lis[-1]
        self.string = string

    def printFASTA(self):
        return '>__XX__' + self.header + '\n' + self.sequence + '\n'


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

    options, args = parser.parse_args()

    global directory
    inputFile = unescape(options.input)
    directory = unescape(options.path) + os.sep + "data"

    os.mkdir(directory)

    if isTabular(inputFile):
        saveMulti(inputFile)
    else:
        saveSingle(inputFile)

    pool = Pool()
    list_of_files = [file for file in os.listdir(directory) if file.lower().endswith(extension)]
    pool.map(prank, list_of_files)
    result = [file for file in os.listdir(directory) if file.lower().endswith(output_extension)]
    with open(directory + os.sep + results, "a") as f:
        for file in result:
            with open(directory + os.sep + file, "r") as r:
                f.write(toData(r.read()) + "\n")

if __name__ == '__main__':
    main()
