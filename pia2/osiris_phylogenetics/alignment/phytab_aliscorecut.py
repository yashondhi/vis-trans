import os
import optparse
import subprocess
from multiprocessing import Pool
import shutil

results_dir = "./data"
results = "results.data"
fasta_extension = ".afa"
alicut_prefix = "ALICUT_"
familyList = []
galaxy_tool_dir = "/home/osiris/bin/"
forbidden_chars = {
    '(': '__rb__',
    ')': '__lb__',
    ':': '__co__',
    ';': '__sc__',
    ',': '__cm__',
    '--': '__dd__',
    '*': '__st__',
    '|': '__pi__',
    ' ': '__sp__'
}


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


def unpackData(families):
    with open(families) as f:
        for line in f:
            seq = Sequence(line)
            with open(results_dir + os.sep + seq.family + fasta_extension, "a") as p:
                p.write(seq.printFASTA())


class Sequence:
    def __init__(self, string):
        lis = string.split('\t')
        self.species = lis[0]
        self.family = lis[1]
        self.name = lis[2]
        self.header = ' '.join(lis[:-1])
        self.sequence = lis[-1]
        self.string = string

    def escapedHeader(self):
        string = self.header
        for key, value in forbidden_chars.iteritems():
            string = string.replace(key, value)
        return string

    def printFASTA(self):
        return '>' + self.escapedHeader() + '\n' + self.sequence + '\n'


def unescapeHeader(header):
    string = header
    for key, value in forbidden_chars.iteritems():
        string = string.replace(value, key)
    return string


def toData(text):
    text = text.split('\n')
    result = ''
    for line in text:
        if '>' in line:
            line = '\n' + unescapeHeader(line.replace('>', "")) + '\t'
        line = line.replace(" ", "\t")
        result += line
    return result[1:]  # Index past the first newline char


def aliscore(input):
    file_name = results_dir + os.sep + input
    # print file_name
    pop = subprocess.Popen(["perl", "-I", galaxy_tool_dir, galaxy_tool_dir + "Aliscore.02.pl", "-i", file_name])
    pop.wait()


def main():
    usage = """%prog [options]
options (listed below) default to 'None' if omitted
    """
    parser = optparse.OptionParser(usage=usage)

    parser.add_option(
        '-i', '--input',
        dest='families',
        action='store',
        type='string',
        metavar="FILE",
        help='Name of input sequences.')

    options, args = parser.parse_args()

    families = unescape(options.families)

    os.mkdir(results_dir)

    unpackData(families)

    list_of_files = [file for file in os.listdir(results_dir) if file.lower().endswith(fasta_extension)]

    pool = Pool()
    pool.map(aliscore, list_of_files)

    alicut = "ALICUT_V2.0_modified.pl"
    shutil.copy(galaxy_tool_dir + alicut, results_dir + os.sep + alicut)
    os.chdir(results_dir)
    pop = subprocess.Popen(["perl", "./" + alicut])
    pop.wait()
    os.chdir("../")

    result = [file for file in os.listdir(results_dir) if file.startswith(alicut_prefix)]
    with open(results_dir + os.sep + results, "a") as f:
        for file in result:
            if file.endswith(fasta_extension):
                with open(results_dir + os.sep + file, "r") as r:
                    f.write(toData(r.read()) + "\n")

if __name__ == '__main__':
    main()
