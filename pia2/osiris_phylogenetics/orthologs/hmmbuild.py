import os
import optparse
import subprocess
from multiprocessing import Pool

directory = "./data"
results = "results.data"
extension = ".afa"
model_extension = ".hmm"
inputFile = ""
index_of_name_in_hmm = 6


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


def toData(text):
    lis = text.split()
    name = lis[index_of_name_in_hmm]
    text = name + "\t" + text.replace("\n", "\\n")
    return text


def hmmbuild(input):
    file_name = directory + os.sep + input
    # print file_name
    # return subprocess.Popen(['hmmbuild', "--informat", "afa", file_name + ".hmm", file_name], stdout=subprocess.PIPE).communicate()[0]  # ./muscle
    pop = subprocess.Popen(['hmmbuild', "--informat", "afa", file_name + ".hmm", file_name])
    pop.wait()


class Sequence:
    def __init__(self, string):
        lis = string.split('\t')
        # print lis
        self.species = lis[0]
        self.family = lis[1]
        self.name = lis[2]
        self.header = ' '.join(lis[:-1])
        self.sequence = lis[-1]
        self.string = string

    def printFASTA(self):
        return '> ' + self.header + '\n' + self.sequence + '\n'


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

    options, args = parser.parse_args()

    global inputFile, directory
    inputFile = unescape(options.input)

    os.mkdir(directory)

    with open(inputFile) as f:
        for line in f:
            seq = Sequence(line)
            with open(directory + os.sep + seq.family + extension, "a") as p:
                p.write(seq.printFASTA())

    pool = Pool()
    list_of_files = [file for file in os.listdir(directory) if file.lower().endswith(extension)]

    pool.map(hmmbuild, list_of_files)

    result = [file for file in os.listdir(directory) if file.lower().endswith(model_extension)]
    with open(directory + os.sep + results, "a") as f:
        for file in result:
            with open(directory + os.sep + file, "r") as r:
                f.write(toData(r.read()) + "\n")

if __name__ == '__main__':
    main()
