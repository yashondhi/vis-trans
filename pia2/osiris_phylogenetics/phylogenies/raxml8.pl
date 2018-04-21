#! /usr/bin/env perl

use strict;
use warnings;
#raxml.pl Galaxy wrapper calls raxml from raxml.xml
#xml file contains:
#raxml.pl [GTR|CAT] [PROT|DNA] [protmodel] [morphmodel] [phylip file] [constraint] [partition] [best_tree?] [invar?] [#bootreps] [outgroup]

##For debugging command line pass, uncomment next
#for (my $i=0; $i < @ARGV; $i++){
#	print "Parameter #$i ".$ARGV[$i]."\n\n";
#}
#exit;

my $rate_het=shift(@ARGV);		#0 rate heterogeneity? value will = GAMMA or CAT
my $datatype = shift(@ARGV);		#1 datatype? True=Protein False=DNA
my $protmodel = shift(@ARGV);		#2 which protein model
my $morphmodel = shift(@ARGV);		#3 which morphology multistate model
my $data_file= shift(@ARGV);		#4 input a phylip file
my $part_file = shift(@ARGV);		#5 optional partition file
my $constraint_tree = shift(@ARGV);	#6 optional constraint tree
my $nboots = shift(@ARGV);		#9 Number of bootstrap reps
my $seed = shift(@ARGV);		#10 Number of bootstrap reps
my $long = shift(@ARGV);		#11 decide whether to do a long or bootstrap call or not, with multiple threads
my $outgroup = shift(@ARGV);		#12 Specify the outgroup
my $model;



# From shell pipeline
#        raxmlHPC-PTHREADS7.2.6 -T $processors -f a -s $data_name.data  -q $data_name.part -m $model -n $data_name -N 100 -x 1234567890 -o Limulus_polyphemus
#        cp RAxML_bestTree.$data_name $data_nameBootBest.tre
#        cp RAxML_bipartitions.$data_name $data_nameBoot.tre

#ADD OPTIONS TO BUILD FULL RAXML COMMANDLINE ARGUMENT

my $build_command;
#First CALL RAXML THROUGH PATH with 8 threads
if($long eq 'Long'){ #Currently both raxml and raxml_long call with 'long'
#	$build_command = "raxmlHPC-PTHREADS-SSE3 -T 8"; This was version 7
	$build_command = "/var/www/galaxy_dev/raxml/standard-RAxML-8.1.1/raxmlHPC-PTHREADS-AVX -T 8";
}else{
	$build_command = "mpirun -np 10 raxmlHPC-MPI-SSE3 ";
}
#Check if bootstrapping is desired
	if($nboots > 0){
		$build_command = $build_command." -f a ";
	}else{
		system "echo '0 bootstraps selected' >> RAxML_bipartitions.galaxy";
		system "echo '0 bootstraps selected' >> RAxML_bipartitionsBranchLabels.galaxy";
		system "echo '0 bootstraps selected' >> RAxML_bootstrap.galaxy";
	}
#Next add call to input phylip file
	$build_command = $build_command." -s ".$data_file;
#Add call to partition file name
	unless($part_file eq 'None'){
		$build_command = $build_command." -q ".$part_file;
	}
#Build substitution model
	if($datatype eq "PROT"){
		$model = "PROT";
	}elsif($datatype eq "DNA"){
		$model = "GTR";
	}elsif($datatype eq "MORPH"){
		$model = "BIN";
	}
	if($rate_het eq "GAMMA"){
		$model = $model."GAMMA";
	}elsif($rate_het eq "CAT"){
		$model = $model."CAT";
	}
	if($datatype eq "PROT"){
		$model = $model.$protmodel;
	}
	$build_command = $build_command." -m ".$model;
#Add multistate morphology model
	unless($morphmodel eq "NULL"){
		$build_command = $build_command." -K ".$morphmodel;
	}
#check constraint tree
	unless($constraint_tree eq "None"){
		$build_command = $build_command." -g ".$constraint_tree;
	}
#IF bootstrap number is zero then skip this section
	if($nboots > 0){
	#N Bootstraps
		$build_command = $build_command." -N ".$nboots;
	#Bootstrap seed
		$build_command = $build_command." -x ".$seed;
	}

#Parsimony seed
	$build_command = $build_command." -p "."1234567";


#name output files galaxy
	$build_command = $build_command." -n galaxy";
#Outgroup
	if(defined $outgroup){
		$build_command = $build_command." -o ".$outgroup;
	}

print "Osiris/Galaxy COMMAND BUILD WAS: $build_command\n";

#Uncomment to actually call raxml
system $build_command;
