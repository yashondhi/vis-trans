#!/usr/bin/perl -w
use strict;

use Bio::SeqIO;


my $datafile = $ARGV[0];
my $outfile = $ARGV[1];

open FILE, ">$outfile" or die "Cannot Write File\n"; 

my $seqio_object = Bio::SeqIO->new(-file => $datafile,'-format' => 'genbank');

while(my $seq_object = $seqio_object->next_seq){
	my $organism = $seq_object->species->binomial();
	$organism =~ s/ /_/g;
	my $accession = $seq_object->id;
	for my $feat_object ($seq_object->get_SeqFeatures) {
	   if ($feat_object->primary_tag eq "CDS") {
	      my $sequence = $feat_object->spliced_seq->seq;
	      	if ($feat_object->has_tag('gene')) {
			for my $name ($feat_object->get_tag_values('product')){
				$name =~ s/ /_/g;
				print FILE $organism."\t".$name."\t".$accession."\t".$sequence."\n";
         		}
		}
	   }elsif ($feat_object->primary_tag eq "misc_RNA") {
	      my $sequence = $feat_object->spliced_seq->seq;
	      	if ($feat_object->has_tag('product')) {
			for my $name ($feat_object->get_tag_values('product')){
				$name =~ s/ /_/g;
				print FILE $organism."\t".$name."\t".$accession."\t".$sequence."\n";
         		}
		}
   	   }

	}
}
close FILE;


