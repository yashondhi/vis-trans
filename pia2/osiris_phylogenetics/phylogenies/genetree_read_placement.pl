#! /usr/bin/env perl

use strict;
use warnings;

##For debugging command line pass, uncomment next
#for (my $i=0; $i < @ARGV; $i++){
#	print "Parameter #$i ".$ARGV[$i]."\n\n";
#}
#exit;

my $newgenes=shift(@ARGV);		#0 new genes to align
my $align = shift(@ARGV);		#1 alignment program to use
my $path = shift(@ARGV);		#2 path to tree and gene data
my $name = shift(@ARGV);		#3 name of gene family

#If $newgenes has not hits, do not do read placement, just write tree with no hits
my $buffer;
my $lines = 0;
open(FILE, $newgenes) or die "Can't open `$newgenes': $!";
while (sysread FILE, $buffer, 4096) {
    $lines += ($buffer =~ tr/\n//);
}
close FILE;

if($lines < 1){
	print "No hits found. Skipping read placement\n Tree copied to output.\n";
	system "cp $path.tre RAxML_labelledTree.EPA_TEST";
	system "cp $path.tre RAxML_originalLabelledTree.EPA_TEST";
}else{

	#First concatenate fasta files and align
	system "cat $newgenes $path.fas > toalign.fas";

	if($align eq "MUSCLE"){
		system "muscle -in toalign.fas -out aligned.fas";
	}
	elsif($align eq "MAFFT") {
        	system "mafft --auto toalign.fas > aligned.fas";
	}
	elsif($align eq "PRANK") {
        	system "prank -d=toalign.fas -o=aligned -f=fasta -F";
		system "mv aligned.2.fas aligned.fas";
	}

	#convert to phylip format, uses seqConverter.pl
	system "perl ../phyloconversion/seqConverterG.pl -daligned.fas -ope -Oaligned.phy";

	system "raxmlHPC-PTHREADS-SSE3 -f v -s aligned.phy -m PROTGAMMAWAG -t $path.tre -n EPA_TEST -T 4";
}

#Now make tab delimited file to use in tab2trees
#open treefile to read tree line
open(TREE, "<","RAxML_labelledTree.EPA_TEST") or die "Can't open RESULT File!";
my $finaltree;
while (<TREE>){
	if($_ =~ /\;/m){
		$finaltree = $_;
		chomp($finaltree);
	}
}	
close TREE;

$name =~ s/ /_/g;
chomp($name);
#remove clade labels
$finaltree =~ s/\[I\d+\]//g;
open(TAB, '>treeout.tab') or die "Can't open File!";         
print TAB $name."\t".$finaltree."\n";
close TAB;
