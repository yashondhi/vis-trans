#!/usr/bin/perl

use strict;

my $fasconcatPath = '/home/galaxy/galaxy_dist/tools/fasconcat/FASconCAT_v1.0.pl';

my $outputFormat = $ARGV[0];
my $limit = $ARGV[1];
my $outFormat;
my @inputFiles;

for(my $i = 2; $i <= $limit; $i++) {
	$inputFiles[$i] = " -f ";
	$inputFiles[$i] = $inputFiles[$i].$ARGV[$i];
}

if($outputFormat == "0") {
	$outFormat = "";
}
elsif($outputFormat == "1") {
	$outFormat = " -p -p";
}
else {
	$outFormat = " -n -n";
}

my $run = qx/$fasconcatPath -s -i "@inputfiles" $outFormat /;

if($outputFormat == "0") {
	qx/cp FcC_smatrix.fas output/;
}
elsif($outputFormat == "1") {
	qx/cp FcC_smatrix.phy output/;
}
else {
	qx/cp FcC_smatrix.nex output/;
}
