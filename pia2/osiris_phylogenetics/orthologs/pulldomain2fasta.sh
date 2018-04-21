#!/bin/sh

hmmsearch -A outfile 
hmmsearch -A outfile $1 $2
esl-reformat fasta outfile >   $3







