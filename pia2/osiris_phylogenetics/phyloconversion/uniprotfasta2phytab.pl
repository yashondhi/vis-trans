#!/usr/bin/perl -w

use strict;

use FindBin;
use lib "$FindBin::Bin/lib";
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Seq;

#inputs
my $infile=shift(@ARGV);
my $partition=shift(@ARGV);
#my $delpipes=shift(@ARGV);
my $species;

my $seqid;
# open infile fasta file
my $in_obj = Bio::SeqIO->new(-file => $infile, '-format' =>'fasta');

while (my $seq = $in_obj->next_seq() ) {
        my $sequence = $seq->seq;
        my @rawid = split(/\|/, $seq->id);
        $seqid = $rawid[1];
#       $seqid = $seq->id;

        $sequence =~ s/\n//g;
        $species = $seq->desc;
        #species Name is after OS=
        $species =~ s/.+OS\=//;
        $species =~ s/.+OS\=//;
        #species Name is before GN= sometimes PE=
        $species =~ s/ GN\=.+//;
        $species =~ s/ PE\=.+//;
        $species =~ s/ /_/g;

        print $species."\t".$partition."\t".$seqid."\t".$sequence."\n";
}
