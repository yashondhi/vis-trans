#!/usr/bin/perl -w
#use strict;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Seq;

my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my $remseqs = $ARGV[2];
my $outfile2 = $ARGV[3];

# open infile fasta file
my $in_obj = Bio::SeqIO->new(-file => $infile, '-format' =>'fasta');
open FILE, ">$outfile" or die $!;
open TOSS, ">$outfile2" or die $!;
my %hash;
my %seqhash;

while (my $seq = $in_obj->next_seq() ) {
	my $sequence = $seq->seq;
	my $seqid = $seq->id;

	my @position = split(/\//, $seqid);
	my $gene = $position[0];
	my @positions = split(/\-/, $position[1]);
	my $from = $positions[0];
	my $to = $positions[1];
	my $domain = [];
	my $domainsorted = [];

	if(($remseqs) && ($gene !~ /seq1\|/)    ){
		print TOSS ">".$gene."\n$sequence\n";
	}else{
		$seqhash{$gene}{$from} = $sequence;

		if(exists $hash{$gene}){
			foreach (@{$hash{$gene}}) {
			    push (@$domain, "$_");
			}
			push (@$domain, $from);
			@$domainsorted = sort { $a <=> $b } @$domain;
			$hash{$gene} = $domainsorted;

		}else{
			push (@$domain, $from);
			$hash{$gene} = $domain;
		}
	}
}


foreach $genename (keys %hash) {
   my $num=0;
   foreach (@{$hash{$genename}}) {
	$num++;
	my $of = @{$hash{$genename}};
	print FILE ">".$genename."_"."$_"."__".$num."of$of\n";
	print FILE $seqhash{$genename}{$_}."\n";
   }
}
