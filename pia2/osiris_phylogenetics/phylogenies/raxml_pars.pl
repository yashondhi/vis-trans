#! /usr/bin/env perl

use strict;
use warnings;
#raxml.pl Galaxy wrapper calls raxml from raxml.xml

#For debugging command line pass, uncomment next 4 lines
#for (my $i=0; $i < @ARGV; $i++){
#	print "Parameter #$i ".$ARGV[$i]."\n\n";
#}
#exit;

my $datatype = shift(@ARGV);		#0 datatype
my $data_file= shift(@ARGV);		#1 input a phylip file
my $part_file = shift(@ARGV);		#2 optional partition file
my $seed = shift(@ARGV);		#3 Number of bootstrap reps
my $outgroup = shift(@ARGV);		#4 Specify the outgroup
my $model;

#ADD OPTIONS TO BUILD FULL RAXML COMMANDLINE ARGUMENT

my $build_command;
#First CALL RAXML THROUGH PATH with 8 threads
	$build_command = "raxmlHPC-PTHREADS-SSE3 ";
#Add Parsimony Option and Thread number for PThreads
	$build_command = $build_command." -y -T 4";
#Next add call to input phylip file
	$build_command = $build_command." -s ".$data_file;
#Add call to partition file name
	unless($part_file eq 'None'){
		$build_command = $build_command." -q ".$part_file;
	}
#model is passed directly with xml
	$model = $datatype;
	$build_command = $build_command." -m ".$model;
#Parsimony seed
	$build_command = $build_command." -p ".$seed;
#name output files galaxy
	$build_command = $build_command." -n parsimony";
#Outgroup
	if(defined $outgroup){
		$build_command = $build_command." -o ".$outgroup;
	}

print "Galaxy COMMAND BUILD WAS: $build_command\n";

#Uncomment to actually call raxml
system $build_command;
