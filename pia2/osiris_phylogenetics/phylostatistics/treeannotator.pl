#!/usr/bin/perl

my $treeannotator_path = '/home/galaxy/pkgs/BEAST172/bin/treeannotator';

my $input = $ARGV[0];
my $burnin = $ARGV[1];
my $Node_heights = $ARGV[2];

my $node_opt;

if($Node_heights eq "0") {
	$node_opt = "keep";
}
elsif($Node_heights eq "1") {
	$node_opt = "median";
}
elsif($Node_heights eq "2") {
	$node_opt = "mean";
}

my $run = qx/$treeannotator_path -heights $node_opt -burnin $burnin $input out.tre 2>log.txt/;

print $run; 
