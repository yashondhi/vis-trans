#!/bin/sh
# set -x
# Uncomment 'set -x' for debug information

#############################################
# 		UCSB HAMSTER - GALAXY HISTORY 		#
# 											#
# Executed when user chooses Galaxy history #
#############################################

# Set your hamster script location here. The directory containing these scripts should be in 
#your path
script="hamstrsearch_local-hmmer3.pl"
# Set your unbuild.py script location here
unbuild="unbuild.py"
# Set your emap2fasta.pl script location here
emap2fasta="emap2fasta.pl"


# 1 - Sequence input file
# 2 - Proteins results output file
# 3 - CDS results output file
# 4 - Screen log
# 5 - Species name
# 6 -  whether to use EST flag D=DNA so use -est flag P=Protein so do not use -est flag in hmmstr call
# 7 - HMM Input from UCSB HMMBUILD
# 8 - MUSCLE data from UCSB MUSCLE
# 9 - Reference Species File
# 10 - Reference Species Name

input=$1
proteins=$2
cdsfile=$3
screenlog=$4
speciesName=$5
datatype=$6
hmm_data=$7
muscle_data=$8
filepath=`pwd`
tail="_prot"
tail2="_temp"

# set flag based on input
if [ $datatype = "P" ];
        then
            	estflag="-protein"
        else
            	estflag="-est"
fi

refspfile=${9}
refsphist=${10}



echo "Protein or EST? : $estflag" >> $screenlog
echo "Reference genome file from galaxy history: $refspfile" >> $screenlog
echo "Reference species genome name: $refsphist" >> $screenlog

# unbuild.py here on $hmm_data
mkdir core
mkdir core/hmm_dir
cp $hmm_data core/core.fa

$unbuild core/hmm_dir core/core.fa
cp core/hmm_dir/hmmlist.txt core/hmmlist.txt

# use formatdb to generate new blastdb from this input file
refsphistGALAXY=$refsphist
mkdir $refsphistGALAXY

cp $muscle_data $refsphistGALAXY/$refsphist$tail2
$emap2fasta $refsphistGALAXY/$refsphist$tail2 $refsphist
cp full.fasta core/core.fa

cp $refspfile $refsphistGALAXY/$refsphist$tail
cd $refsphistGALAXY
formatdb -t $refsphist -i $refsphist$tail -n $refsphist$tail

echo "*** Direcotry Structure of Ref. Genome ***" >> $screenlog
ls -l >> $screenlog
echo >> $screenlog

cd $filepath

# script execution
$script -sequence_file=$1 $estflag -taxon=$5 -hmmset=core -refspec=$refsphistGALAXY -galaxyout=$2 -2galaxyout=$cdsfile 2>log.txt >> $screenlog
