#!/usr/bin/perl -w

use strict;

use FindBin;
use lib "$FindBin::Bin/lib";
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Seq;

#inputs
my $infile=shift(@ARGV);
my $species=shift(@ARGV);
my $partition=shift(@ARGV);
my $delpipes=shift(@ARGV);
my $fromfasta;
#for debugging xml input
#print "$infile $species $partition $delpipes\n";
#exit;

if($species eq "from fasta"){
	$fromfasta=1;
}
my $seqid;
# open infile fasta file
my $in_obj = Bio::SeqIO->new(-file => $infile, '-format' =>'fasta');

#no warnings 'uninitialized';	#Was getting error on one fasta for uninitialized sequences. Never could track down why and used this as a workaround
while (my $seq = $in_obj->next_seq() ) {
	my $sequence = $seq->seq;
	$seqid = $seq->id;
	if($delpipes eq 'yes'){
		$seqid =~ s/\|/_/g;
	}
	if($fromfasta){
		$species = $seqid;
	} 
	$sequence =~ s/\n//g;
	print $species."\t".$partition."\t".$seqid."\t".$sequence."\n";
}
