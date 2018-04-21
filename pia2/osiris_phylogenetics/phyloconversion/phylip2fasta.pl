#!/usr/bin/perl -w


my $infile = $ARGV[0];
my $interlv = $ARGV[1];
my $idlen = $ARGV[2];
my $outfile = $ARGV[3];

#open (OUT, ">$outfile");

use strict;
    use Bio::AlignIO;
    use Bio::SimpleAlign;
    #you can set the name length to something other than the default 10
    #if you use a version of phylip (hacked) that accepts ids > 10
    my $outstream = Bio::AlignIO->new(-format  => 'fasta',
                                        -fh      => \*STDOUT );

    # convert data from one format to another
    my $phylipstream     =  Bio::AlignIO->new(-interleaved => $interlv,
					  -format => 'phylip',
					  -file => '<'.$infile,
					  -idlength=>$idlen );
  while( my $aln = $phylipstream->next_aln ) {
        $outstream->write_aln($aln);
    }

