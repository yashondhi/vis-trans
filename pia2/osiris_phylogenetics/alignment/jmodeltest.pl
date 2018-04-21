#!/usr/bin/perl

use strict;
use warnings;
use Cwd;

my $jmodeltest_path = '/home/galaxy2/pkgs/jmodeltest'; # where the JAR package is located
my $jmodeltest_tool = '/home/galaxy2/galaxy-dist/tools/osiris_phylogenetics/alignment'; # where this tool is installed for galaxy

my $galaxyPath = getcwd();

my $input = $ARGV[0];
my $likelihoodStyle = $ARGV[1];
my $likelihood = $ARGV[2];
my $criterion = $ARGV[3];
my $extension = $ARGV[4];

# open increment file
open my $file, '<', $jmodeltest_tool."\/increment.txt";
	my $increment = <$file>;
	$increment = int($increment);
close $file;

# get the current increment
my $temp = $increment;

# update the increment
open(UPDATE, '>'.$jmodeltest_tool."\/increment.txt");
	$increment = $increment + 1;
	print UPDATE $increment;
close(UPDATE);

chdir("$jmodeltest_path");

if($likelihoodStyle eq "-t") {
	#only need to copy input file
	qx/cp $input input.$temp.$extension/;
	
	#print qx/ls/;
	
	qx/java -jar jModelTest.jar -d input.$temp.$extension $likelihoodStyle $likelihood -$criterion > $galaxyPath\/output.txt 2> $galaxyPath\/err_log.txt/;		
	qx/rm input.$temp.*/;
}
elsif($likelihoodStyle eq "-u") {
	#copy input file
	qx/cp $input input.$temp.$extension/;
	#copy likelihood tree
	qx/cp $likelihood likelihood.$temp.tre/;
	qx/java -jar jModelTest.jar -d input.$temp.$extension $likelihoodStyle likelihood.$temp.tre -$criterion > $galaxyPath\/output.txt 2> $galaxyPath\/err_log.txt/;		

	# clean up
	qx/rm input.$temp.*/;
	qx/rm likelihood.$temp.tre/;

}

chdir("$galaxyPath");
