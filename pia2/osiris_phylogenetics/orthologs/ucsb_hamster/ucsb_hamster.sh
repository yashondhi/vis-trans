#!/bin/sh

#############################################
# 				UCSB HAMSTER				#
# 											#
# Executes Hamster with given XML parameter #
#############################################

#the name of the script is here
#script="/home/osiris/galaxy-dist/tools/osiris/orthologs/ucsb_hamster/hamstrsearch_local-hmmer3.pl"
script="hamstrsearch_local-hmmer3.pl"
 


#Variables input from xml
# 1 - Sequence input file
# 2 - Proteins Results Output file
# 3 - cds Results Output file
# 4 - Screen Log
# 5 - Species Name
# 6 - Whether to use -est flag (if D) or not (if P)
# 7 - Core ortholog name
# 8 - Base path for local core orthologs
# 9 - Base path for local reference blast database
#10 - Reference genome

input=$1
proteins=$2
cdsfile=$3
screenlog=$4
species=$5
datatype=$6
core=$7
corepath=$8
blastpath=$9
genome=${10}

echo "ucsb_hamster.sh script parameters" >> $screenlog
echo "Core ortholog name is $core " >> $screenlog
echo "Reference genome name is $genome " >> $screenlog
echo "Species name is $species " >> $screenlog
echo "Datatype $datatype " >> $screenlog

#set flag based on input
if [ $datatype = "P" ];
        then
            	estflag="-protein"
        else
            	estflag="-est"
        fi

# First copy hmm's to working directory
# Currently copies from Data directory
mkdir $core
mkdir $core/hmm_dir
cp -r $corepath/* ./$core/
echo "cp $corepath/hmm_dir/* ./$core/hmm_dir/" >> $screenlog
cp -r $corepath/hmm_dir/* ./$core/hmm_dir/

# Currently copies from data directory
mkdir $genome
cp $blastpath/* ./$genome/

# Now call the actual Hamster Script
#$script -sequence_file=$input $estflag -taxon=$species -hmmset=$core -refspec=$genome -galaxyout=$proteins -2galaxyout=$cdsfile >> $screenlog 2>log.txt
$script -sequence_file=$input $estflag -taxon=$species -hmmset=$core -refspec=$genome -galaxyout=$proteins -2galaxyout=$cdsfile >> $screenlog

