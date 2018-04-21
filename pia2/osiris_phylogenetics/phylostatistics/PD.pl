#!/usr/bin/perl -w

use strict;

#use FindBin;
#use lib "$FindBin::Bin/lib";
use Bio::TreeIO;
use Bio::Tree::Tree;

###this script will find the phylogenetic distance between two species
#input is a tree, output filename, and table with pairwise distances
#usage:
#PD.pl <pairsTable> <treefile> <outfile> <yes|no>
# parse in newick/new hampshire format
my @species1;
my @species2;


my $half=$ARGV[3];
my $divtimebool;
if($half eq 'yes'){
	$divtimebool=1;
}elsif($half eq 'no'){
	$divtimebool=0;
}else{
	die "Argument must contain yes or no for divergence times\n";
}
my $outfile = $ARGV[2];
open(OUT, ">$outfile") or die("Couldn't open output file $ARGV[2]\n");


my $pairsfile = $ARGV[0];
open(PAIRS, "$pairsfile") or die("Couldn't open input file $ARGV[0]\n");
while (<PAIRS>) {
        chomp;
        my $sp1;
        my $sp2;
        ($sp1, $sp2) = split("\t");
        push(@species1, $sp1);
        push(@species2, $sp2);
}

my $treefile = $ARGV[1];

for(my $i=0; $i < @species1; $i++){
        print OUT $species1[$i]."\t".$species2[$i];
        open(TREE, "$treefile") or die("Couldn't open output file $ARGV[1]\n");

        my $treeio = new Bio::TreeIO('-format' => 'newick',
                                   '-file'   => $treefile);

        while(my $tree = $treeio->next_tree){;
                my $node1 = $tree->find_node(-id => $species1[$i]);
                my $node2 = $tree->find_node(-id => $species2[$i]);
                my $distances = $tree->distance(-nodes => [$node1,$node2]);

                #ADD OPTION FOR DIVIDING BY 2 FOR DIVERGENCE TIMES
		if($divtimebool==1){
                	$distances = $distances/2 ;
		}
                print OUT "\t".$distances;
        }
print OUT "\n";
close(TREE);
}

close(PAIRS);
close(OUT);
