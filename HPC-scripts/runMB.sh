#!/bin/bash

#from commmand line: bsub -n 8 -q PQ_liberles /home/liberles/bin/runMB.sh nexusfile $PWD
module load mrbayes/3.2.2

#list of files (complete names ,with extension)
set nexusfile = $1 # nexus file for Mr. Bayes run
set workdir = $2 #PWD

cd $workdir/
mpirun -np 8 mb $nexusfile >& ${nexusfile}_MB.log
