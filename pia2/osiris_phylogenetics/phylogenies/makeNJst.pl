#!/usr/bin/env perl

#This script generates an R script to call NJst
#input is a table with treename<tab>newick tree
use strict;
use Bio::TreeIO;

my $filename = $ARGV[0];
my $outfile = $ARGV[1];
open FILE, $filename or die $!;


my @splitline;

print "require(phybase);\n";
print "genetrees<-c(";
my $counter=0;
my $tree;
while (<FILE>) {
        chomp;
        #get a line from the data file
        my $currentinput = "$_";
	@splitline = split(/\t/);
	my $treename= $splitline[0];
	$tree = $splitline[1];
	unless($counter==0){
		print ", ";
	}
	$counter++;
        print "'$tree'";
}
print ")\n"; #close genetree vector
print "taxaname<-c(";
my $spnum = tree2spList($tree);
print ")\nspname<-taxaname\n";
print "species.structure<-matrix(0,$spnum,$spnum)\n";
print "diag(species.structure)<-1\n";
print "\n";
print "result<-NJst(genetrees,taxaname,spname,species.structure)\n";
print "write(result, file='$outfile')\n";
close FILE;





#This script requires phybase R package
#NJst is a function used as follows
#	genetrees<-c("(A:0.004,(B:0.003,(C:0.002,(D:0.001,E:0.001)
#		:0.001):0.001):0.001);","(A:0.004,(B:0.003,(E:0.002,(D:0.001,C:0.001):0.001):0.001):0.001);","(A:0.004,(B:0.003,(C:0.002,(D:0.001,E:0.001):0.001):0.001):0.001);")
#     taxaname<-c("A","B","C","D","E")
#     spname<-taxaname
#     species.structure<-matrix(0, 5, 5)
#     diag(species.structure)<-1
#     
#     NJst(genetrees,taxaname, spname, species.structure)



sub tree2spList {
	my $treefile=shift;

	my ($charactername, $characterstate); 
	my ($call, $sp_id, $char_id);

	#Open treefile and get taxon names from tree
	my $stringfh;
	open($stringfh, "<", \$treefile);

	my $input = Bio::TreeIO->new(-format => 'newick', -fh => $stringfh); 
	my $tree = $input->next_tree; 

	my @taxa = $tree->get_leaf_nodes; 
	my @names = map { $_->id } @taxa;

	my $count=0;
	foreach(@names){
		my $treespecies = $_;
		$treespecies =~ s/^\s+|\s+$//g ;	#Trim leading and trailing whitespace
		unless($count==0){
			print ",";
		}
		print "'$treespecies'";
		$count++
	}
	return $count;
}	#end of tree2spList subroutine
