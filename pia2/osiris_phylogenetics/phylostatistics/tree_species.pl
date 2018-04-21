#!/usr/bin/perl
use strict;

my $infile=$ARGV[0];
my $treefile=$ARGV[1];

open IN, $infile or die "Cannot open $infile\n";
open TREE, $treefile or die "Cannot open $treefile\n";

my $tree = <TREE>;
close(TREE);

while(<IN>){
	my $curspecies = $_;
	chomp($curspecies);
	$curspecies =~ s/ /\_/g ;
	if($tree =~ m/$curspecies/){
		#match
	}else{
		print $curspecies."\n";
	}	
}
