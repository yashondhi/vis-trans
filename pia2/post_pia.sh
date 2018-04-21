#!/bin/bash
### This collection of scripts is written in bash, perl, and python.
###
### Make sure dependencies are installed and configured before running.
### Dependencies include: BioPerl, NumPy, SciPy, statsmodels (for python).
###
### With the exception of cleanhits.py and the customized long_branch_finder2.py, 
### all other scripts can be found as part of the osiris_phylogenetics galaxy module
### avalaible here: https://bitbucket.org/osiris_phylogenetics/osiris_phylogenetics
### 
### Special thanks to Todd Oakley's lab at UCSB for the original galaxy scripts on 
### which we have based this pipeline.
###
###  Jorge L. Perez-Moreno, Danielle DeLeo
###  CRUSTOMICS Lab at Florida International University 04/05/2017



###  Identify tips longer than ___ median absolute deviations of the tree's branch lengths. Modify number to change MAD multiplier.
###  long_branch_finder2.py can be reverted to calculate standard deviations instead of MADs. More info as comments in the script.

#python ~/apps/osiris_phylogenetics/phylogenies/long_branch_finder2.py treeout.tab 4 > hits_to_prune.lit
long_branch_finder2.py treeout.tab 4 > hits_to_prune.list


### Clean the output of the long branch finder to avoid conflicts downstream.

#python ~/apps/pia/cleanhits.py hits_to_prune.list > hits_to_prune.clean.list
cleanhits.py hits_to_prune.list > hits_to_prune.clean.list


### Fix PIA's allhits.tab into proper phytab, then removes the | from old assemblies.

awk -F '\t' '{print $1"\t"$3"\t"$2}' allhits.tab > allhits.fixed.tab
sed -ie "s/|/_/g" allhits.fixed.tab


### Remove entries from PIA results phytab file that match to a list.

#python ~/apps/osiris_phylogenetics/phyloconversion/prune_phytab_using_list.py allhits.fixed.tab hits_to_prune.clean.list discard > allhits.pruned.tab
prune_phytab_using_list.py allhits.fixed.tab hits_to_prune.clean.list discard > allhits.pruned.tab


### Convert back to FASTA


awk '{print ">"$1"_"$2"\n"$3}' allhits.pruned.tab > allhits.pruned.fasta


### Remove duplicated sequences resulting from translation of similar isoforms.

usearch -cluster_fast allhits.pruned.fasta -sort length -id 1.00 -threads 8 -centroids PIA.results.fasta
