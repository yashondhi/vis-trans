#!/usr/bin/perl -w

use strict;

use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Seq;

#inputs
my $infile=shift(@ARGV);
my $outfile=shift(@ARGV);
my $deloutfile=shift(@ARGV);
my $percent=shift(@ARGV);

my $seqid;
my $newnumbers=1; 		#for sequential renumbering of header

open FILE, ">$outfile" or die $!;
open DELFILE, ">$deloutfile" or die $!;


# open infile fasta file to get average length
my $in_obj = Bio::SeqIO->new(-file => $infile, '-format' =>'fasta');

my $seqcount;
my $seqsum = 0;
my $avelen;
while (my $seq = $in_obj->next_seq() ) {
	my $sequence = $seq->seq;
	if($sequence){
		my $seqlen = length($sequence);
		$seqcount++;
		$seqsum = $seqsum + $seqlen;
	}
}
$avelen = $seqsum/$seqcount;
print "AVE= $avelen \n";


# open infile fasta file to get average length
$in_obj = Bio::SeqIO->new(-file => $infile, '-format' =>'fasta');

while (my $seq = $in_obj->next_seq() ) {
	my $sequence = $seq->seq;
	$seqid = $seq->id;
	if($sequence){
		$sequence =~ s/\n//g;
		$sequence =~ tr/a-z/A-Z/;
		my $seqlen = length($sequence);
	
		if($seqlen > ($avelen * ($percent/100) ) ){
			print FILE ">";
			print FILE $seqid." ".$seq->desc."\n".$sequence."\n";
		}else{
			#print "Writing sequence of $seqlen to DELFILE\n";
			print DELFILE ">";
			print DELFILE $seqid." ".$seq->desc."\n".$sequence."\n";
		}		
	}
}

close FILE;
close DELFILE;
