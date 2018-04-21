#! /usr/bin/perl -w

use strict;
use warnings;
use Bio::Seq;
use Bio::AlignIO;
use Bio::Tools::IUPAC;

#Written by Todd H. Oakley UCSB

#Obtain Arguments
my $infile1=shift(@ARGV);		#0 input file
my $threshold=shift(@ARGV);		#6 data outfile

my %HoGF;
my $line = "";
my @specieslist;


#INPUT ALL DATA INTO HASH
open(INFILE1, $infile1);
foreach $line(<INFILE1>) {
	chop($line);
	my $getline = $line;
	my @column = split(/\t/, $getline);
	my $species = $column[0];
	my $genefamily = $column[1];
	my $genename = $column[2];
	my $sequence = $column[3];

	$HoGF{$genefamily}{$species}=$sequence;
	push(@specieslist,$species);
}
@specieslist = uniq(@specieslist); 	#uniq subroutine defined at end

foreach my $gfkey (sort keys %HoGF) {
	my $fastalikestring="";
	for(my $j=0; $j<@specieslist; $j++){
	               if(exists $HoGF{$gfkey}{$specieslist[$j]}){
			$fastalikestring = $fastalikestring.">".$specieslist[$j]."_".$gfkey."\n".$HoGF{$gfkey}{$specieslist[$j]}."\n";
		}
	}
	print $gfkey."\t".cons($fastalikestring, $threshold)."\n";
}
exit;

sub uniq {
    return keys %{{ map { $_ => 1 } @_ }};
}

sub cons {
	my ($fastalikestring, $threshold) = @_;

	my $alnio = Bio::AlignIO->new(-string => $fastalikestring, -format => 'fasta');

	while (my $aln = $alnio->next_aln) {
	    my $ci = $aln->consensus_iupac;
	    my $cs = $aln->consensus_string($threshold);
	    return($cs);
	}
}
