#!/bin/bash

# Helper script to change output file names then run BEAST.

# Usage: ./beast.sh STDOUT_LOG XMLCONFIG BEAGLESSE

newxml=$(python26 /home/galaxy/galaxy-dist/tools/osiris/phylogenies/beast.py ${2} $(pwd))
java -jar -Xms4096m -Xmx8192m /home/galaxy/pkgs/BEAST172/lib/beast.jar -overwrite ${3} -threads 8 ${newxml} > ${1} 2>&1
