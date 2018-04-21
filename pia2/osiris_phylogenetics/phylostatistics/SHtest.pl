#! /usr/bin/perl -w

use strict;
use warnings;

##For debugging command line pass, uncomment next
#for (my $i=0; $i < @ARGV; $i++){
#	print "Parameter #$i ".$ARGV[$i]."\n\n";
#}
#exit;

my $datatype = shift(@ARGV);		#0 datatype
my $data_file= shift(@ARGV);		#1 input a phylip file
my $part_file = shift(@ARGV);		#2 optional partition file
my $best_tree = shift(@ARGV);		#3 best tree for SH comparison
my $alt_trees = shift(@ARGV);		#4 Alternative tree(s) for SH comparison
my $model;

#ADD OPTIONS TO BUILD FULL RAXML COMMANDLINE ARGUMENT

my $build_command;
#First CALL RAXML THROUGH PATH with 8 threads
	$build_command = "raxmlHPC-PTHREADS-SSE3 ";
#Add SH Test Option and Thread number for PThreads
	$build_command = $build_command."-f h -T 4";
#Next add call to input phylip file
	$build_command = $build_command." -s ".$data_file;
#model is passed directly with xml
	$model = $datatype;
	$build_command = $build_command." -m ".$model;
#Add call to partition file name
	unless($part_file eq 'None'){
		$build_command = $build_command." -q ".$part_file;
	}
#Next add call to input best tree file
	$build_command = $build_command." -t ".$best_tree;
#Next add call to input best tree file
	$build_command = $build_command." -z ".$alt_trees;
#name output files galaxy
	$build_command = $build_command." -n SH";

print "Galaxy COMMAND BUILD WAS: $build_command\n";

#Uncomment to actually call raxml
system $build_command;

