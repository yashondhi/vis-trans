#! /usr/bin/perl -w

use strict;
use warnings;
#gblocks.pl [fasta file]

my $infile=shift(@ARGV);
my $datatype=shift(@ARGV);
my $gaps=shift(@ARGV);
my $size=shift(@ARGV);
my $outfileloc=shift(@ARGV);
my $htmlfileloc=shift(@ARGV);




##For debugging command line pass, uncomment next
#for (my $i=0; $i < @ARGV; $i++){
#	print "Parameter #$i ".$ARGV[$i]."\n\n";
#}

system "Gblocks $infile $datatype $gaps -b4=$size";

#Gblocks requires output from $input.fas to be written to $input.fas-gb
#Copy that file to gout where galaxy expects to find the output
my $outfile = $infile."-gb";
my $htmlfile = $outfile.".htm";
system "cat $outfile > $outfileloc";
system "cat $htmlfile > $htmlfileloc";
exit;
