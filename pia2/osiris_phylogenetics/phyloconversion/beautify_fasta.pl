#!/usr/bin/perl -w

use strict;

use FindBin;
use lib "$FindBin::Bin/lib";
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Seq;

#inputs
my $infile=shift(@ARGV);
my $outfile=shift(@ARGV);
my $delpipes=shift(@ARGV);
my $convgi=shift(@ARGV);
my $delslash=shift(@ARGV);
my $renumber=shift(@ARGV);
my $space=shift(@ARGV);
my $truncate=shift(@ARGV);
my $suffix=shift(@ARGV);

my $seqid;
my $newnumbers=1; 		#for sequential renumbering of header
# open infile fasta file
my $in_obj = Bio::SeqIO->new(-file => $infile, '-format' =>'fasta');
open FILE, ">$outfile" or die $!;

while (my $seq = $in_obj->next_seq() ) {
	my $sequence = $seq->seq;
	$seqid = $seq->id;

	if($delslash eq 'yes'){
		$seqid =~ s/\\/_/g;
	}
	if($convgi eq 'yes'){
		$seqid =~ s/gi\|/gi_/g;
	}
	if($delpipes eq 'yes'){
		$seqid =~ s/\|/ /g;
	}
	if($suffix eq 'none'){
		#do nothing
	}elsif($suffix eq '.1'){
		$seqid = $seqid.".1";
	}elsif($suffix eq '.2'){
		$seqid = $seqid.".2";
	}
	$sequence =~ s/\n//g;
	$sequence =~ tr/a-z/A-Z/;
	print FILE ">";
	if($renumber eq 'yes'){
		print FILE $newnumbers;
		if($space eq 'yes'){
			print FILE " ";
		}
		$newnumbers++;
	}
	print FILE $seqid;
		if($truncate eq 'no'){
			print FILE " ".$seq->desc;
		}
	print FILE "\n".$sequence."\n";
}
