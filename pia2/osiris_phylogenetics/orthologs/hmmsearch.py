import os
import optparse
import subprocess
from multiprocessing import Pool

results_dir = "./data"
results = "results.data"
result_extension = ".out"
model_extension = ".hmm"
database = ""


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


def unpackData(models):
    with open(models) as f:
        for line in f:
            hmm = HMM(line)
            with open(results_dir + os.sep + hmm.name + model_extension, "a") as p:
                # print(hmm.model)
                p.write(hmm.model)


class HMM:
    def __init__(self, string):
        lis = string.split('\t')
        # print lis
        self.model = self.restoreNewLines(lis[1])
        self.name = lis[0]

    def restoreNewLines(self, string):
        return string.replace('\\n', '\n')


def toData(text):
    # lis = text.split()
    # name = lis[index_of_name_in_hmm]
    # text = name + "\t" + text.replace("\n", "\\n")
    # text = text.replace("\n", "\\n")
    return text


def hmmsearch(input):
    file_name = results_dir + os.sep + input
    # print file_name
    # return subprocess.Popen(['hmmbuild', "--informat", "afa", file_name + ".hmm", file_name], stdout=subprocess.PIPE).communicate()[0]  # ./muscle
    pop = subprocess.Popen(['hmmsearch', "-o", file_name + result_extension, file_name, database])
    pop.wait()


def main():
    usage = """%prog [options]
options (listed below) default to 'None' if omitted
    """
    parser = optparse.OptionParser(usage=usage)

    parser.add_option(
        '-i', '--hmm',
        dest='hmm',
        action='store',
        type='string',
        metavar="FILE",
        help='Name of input hmm models.')

    parser.add_option(
        '-d', '--database',
        dest='database',
        action='store',
        type='string',
        metavar="FILE",
        help='Name of sequence database.')

    options, args = parser.parse_args()

    global database
    models = unescape(options.hmm)
    database = unescape(options.database)

    os.mkdir(results_dir)

    unpackData(models)

    list_of_files = [file for file in os.listdir(results_dir) if file.lower().endswith(model_extension)]

    pool = Pool()
    pool.map(hmmsearch, list_of_files)

    result = [file for file in os.listdir(results_dir) if file.lower().endswith(result_extension)]
    with open(results_dir + os.sep + results, "a") as f:
        for file in result:
            with open(results_dir + os.sep + file, "r") as r:
                f.write(toData(r.read()) + "\n")

if __name__ == '__main__':
    main()
