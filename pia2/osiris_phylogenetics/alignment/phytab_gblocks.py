#!/usr/bin/python2.7

import phytab_splitter
import subprocess
import argparse
import zipfile
import os

parser = argparse.ArgumentParser(description='Run Gblocks on Phytab')
parser.add_argument('-i', '--input', required=True, dest='input_filename')
parser.add_argument('-t', '--type', required=True, dest='input_type')
parser.add_argument('-b5', '--gap-allowance', required=True, dest='gap_allowance')
parser.add_argument('-b4', '--minblock', type=int, required=True, dest='block_size', help='Minimum Length Of A Block')
parser.add_argument('-b3', '--maxcontigpos', type=int, required=True, dest='max_pos', help='Maximum Number Of Contiguous Nonconserved Positions')
parser.add_argument('-o', '--output', required=True, dest='output_filename')
parser.add_argument('--phyout', required=True, dest='phytab_output_filename')
args = parser.parse_args()

zip = zipfile.ZipFile(args.output_filename, 'w')

if args.input_type == 'd':
    sp = phytab_splitter.Splitter(dnaFile=args.input_filename, minSeqSize=2)
    sp.generateFasta()
    info = sp.getDNAInfo()
else:
    sp = phytab_splitter.Splitter(protFile=args.input_filename, minSeqSize=2)
    sp.generateFasta()
    info = sp.getProteinInfo()

output_info = info

for key in info:
    if args.input_type == 'd':
        print('processing ' + key + 'DNA.fas')
        p = subprocess.Popen(['Gblocks', key + 'DNA.fas', '-t=d', '-b5=' + args.gap_allowance, '-b4',
                              str(args.block_size), '-b3', str(args.max_pos)])
        p.communicate()
        zip.write(key + 'DNA.fas-gb')
        zip.write(key + 'DNA.fas-gb.htm')
        with open(key + 'DNA.fas-gb', 'r') as f:
            lines = f.readlines()
            curr_sequence = []
            species_count = 0
            line_count = 0 # first line is line 0
            for line in lines:
                line_count += 1
                if (line[0] != '>'):
                    curr_sequence.append(line)
                if line_count == (len(lines) - 1):
                    curr_sequence.append(lines[line_count])
                    output_info[key][species_count][2] = ''.join(''.join(curr_sequence).split())
                    break
                if (lines[line_count][0] == '>'):
                    output_info[key][species_count][2] = ''.join(''.join(curr_sequence).split())
                    curr_sequence = []
                    species_count += 1
    else:
        p = subprocess.Popen(['Gblocks', key + 'Prot.fas', '-t=d', '-b5=' + args.gap_allowance, '-b4',
                              str(args.block_size), '-b3', str(args.max_pos)])
        p.communicate()
        try:
            zip.write(key + 'Prot.fas-gb')
            zip.write(key + 'Prot.fas-gb.htm')
            with open(key + 'Prot.fas-gb', 'r') as f:
                lines = f.readlines()
                curr_sequence = []
                species_count = 0
                line_count = 0
                for line in lines:
                    line_count += 1
                    if (line[0] != '>'):
                        curr_sequence.append(line)
                    if line_count == (len(lines) - 1):
                        curr_sequence.append(lines[line_count])
                        output_info[key][species_count][2] = ''.join(''.join(curr_sequence).split()) 
                        break
                    if (lines[line_count][0] == '>'):
                        output_info[key][species_count][2] = ''.join(''.join(curr_sequence).split())
                        curr_sequence = []
                        species_count += 1
        except OSError as e:
                print('I/O error({0}): {1}. Most likely Gblocks was unable to produce the file {2}. Ensure your data is consistent '
                      'and in the right format and try again.'.format(e.errno, e.strerror, key + 'Prot.fas-gb'))
                continue


with open(args.phytab_output_filename, 'w') as f_out:
    out_lines = []
    out_line = ''
    for key in output_info:
        for i in range(0, len(output_info[key])):
            out_line = output_info[key][i][0] + '\t' + key + '\t' + output_info[key][i][1]  + '\t' + \
                       output_info[key][i][2] + '\n'
            f_out.write(out_line)

