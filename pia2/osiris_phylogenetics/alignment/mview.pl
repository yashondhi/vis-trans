#!/usr/bin/perl

my $input = $ARGV[0];
my $dna = $ARGV[1];

if ($dna eq 'dna'){
	$dna = '-DNA';
}else{
	$dna = '';
}
my $run = qx/mview -in pearson $dna -bold -coloring group -html head $input/;
print $run;
