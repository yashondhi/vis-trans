#!/usr/bin/env perl
use strict;

#use FindBin;
#use lib "$FindBin::Bin/lib";
use Bio::DB::GenBank;
use Bio::SeqIO;

my $accession = $ARGV[0];
my $datatype = $ARGV[1];
my $outtype = $ARGV[2];
my $outfile = $ARGV[3];


        my $qry_string .= $accession."[accession]";
        my $fh = Bio::SeqIO->newFh(-format=>$outtype, -file=>">$outfile");

        my $GBseq;
        my $gb = new Bio::DB::GenBank;
        my $query = Bio::DB::Query::GenBank->new
                (-query   =>$qry_string,
                 -db      =>$datatype);

        my $count;
        my $species;
        my $seqio = $gb->get_Stream_by_query($query);
        while( defined ($GBseq = $seqio->next_seq )) {
                my $sequence = $GBseq;   # read a sequence object
                print $fh $sequence; # write a sequence object
        }

exit;
