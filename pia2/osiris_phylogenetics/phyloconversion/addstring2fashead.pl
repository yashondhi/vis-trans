#!/usr/bin/perl -w

use strict;

use FindBin;
use lib "$FindBin::Bin/lib";
#use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Seq;

# open infile fasta file
  my $in_obj      = Bio::SeqIO->new(-file => $ARGV[0], '-format' =>'fasta');

  my $currentinput =  $ARGV[1];

# grab sequence object
        while (my $seq = $in_obj->next_seq()  ) {
        	my $seq_obj = $in_obj;
		print ">".$currentinput.$seq->id."\n";
		print $seq->seq."\n";
        }
