#!/bin/bash

mpirun -np 16 raxmlHPC-MPI -s ${1} -b 12345 -N 50 -m PROTGAMMAJTT -n raxtest -p 12345
