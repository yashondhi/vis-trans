#!/usr/bin/env python

"""
Trims bootstrap probabilities below a user-defined threshold
"""

import argparse


def trim(input_file, probability_ceiling):
    lines_to_save = []
    with input_file as f_in:
        for line in f_in.readlines():
            if line.strip() != '':
                probability = line.split()[-1]
                if float(probability) >= probability_ceiling:
                    lines_to_save.append(line)
        f_in.close()
    return lines_to_save

arg_parser = argparse.ArgumentParser(description='Trims bootstrap probabilities'
                                                 ' below a threshold value.')
arg_parser.add_argument('--input', '-i', metavar='FILENAME', dest='in_file',
                        type=argparse.FileType('r'), required=True,
                        help='Location of input file; format should be tabular '
                             ' with probabilities in decimal form listed in the '
                             'last column')
arg_parser.add_argument('--output', '-o', metavar='FILENAME', dest='out_file',
                        type=argparse.FileType('w'), required=True,
                        help='Location of output file, output file like input, '
                             'except lines with probabilities below threshold '
                             'are removed')
arg_parser.add_argument('--prob_limit', '--prob', '-p', metavar='PROBABILITY',
                        dest='prob_limit', type=float, required=True,
                        help='Decimal probability value; lines in input file '
                             'that have probabilities below this number will '
                             'not be written to the output file')
args = arg_parser.parse_args()

out_lines = trim(args.in_file, args.prob_limit)
with args.out_file as f:
    for item in out_lines:
        f.write(item)
    f.close()
