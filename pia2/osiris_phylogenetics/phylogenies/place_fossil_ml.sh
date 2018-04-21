#!/bin/bash
datafile=$1	#binary datafile
model=BINGAMMA
morphmodel=MK
tree=$2
reps=100

raxmlHPC-PTHREADS-SSE3 -T 8 -f u -s $datafile -K $morphmodel -m $model -n galaxy -t $tree -N $reps
raxmlHPC-PTHREADS-SSE3 -T 8 -f v -s $datafile -K $morphmodel -m $model -a RAxML_weights.galaxy -n fossil_weights -t $tree -N $reps
