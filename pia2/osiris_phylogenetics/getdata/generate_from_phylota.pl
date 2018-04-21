#!/usr/bin/env perl
use strict;
use LWP::Simple;
use Bio::SeqIO;

my $infile = $ARGV[0];
my $database = $ARGV[1];
my $outfile = $ARGV[2];
my $treefile = $ARGV[3];
my $phytabfile = $ARGV[4];

open(IN, "$infile") or exit;
open(OUT, ">$outfile") or exit;
open(TREES, ">$treefile") or exit;
open(DATA, "<$database") or exit;
my @data=<DATA>;
close(DATA);
my @foundtrees;
my @alltrees;
while(<IN>){
        my $line = $_;
	chomp($line);
	$line =~ s/ /_/g;
	print "Finding trees with $line ....";
	@foundtrees = grep /$line/, @data;
	my $numtrees = scalar @foundtrees;
	print " $numtrees Trees\n";
	if($numtrees == 0){
		my @genus = split(/_/,$line);
		@foundtrees = grep /$genus[0]/, @data;
		print "\tTrying genus $genus[0]";
		my $numtrees = scalar @foundtrees;
		print " $numtrees genus Trees\n";
	}
	push(@alltrees,@foundtrees);
}

@alltrees = uniq(@alltrees);
print TREES @alltrees;

#get fasta files for trees
for(my $i=0;$i < @alltrees; $i++){
	my @tablines = split(/\t/,$alltrees[$i]);
	my @tici = split(/_/, $tablines[0]);
	my $ti = $tici[0];
	my $ci = $tici[1];
	my $addstring = $ti.$ci."_";
	$ti =~ s/ti//;
	$ci =~ s/ci//;
	my $fastafile = getfastafromphylota($ci,$ti);
	#Add TI_CI_ to each fastaheader
	$fastafile =~ s/\>/\>$addstring/g;
	print OUT $fastafile;
}
close(IN);
close(OUT);
close(TREES);

#Now convert fasta file to phytab file and write
open(PHYTAB, ">$phytabfile") or exit;
# open infile fasta file
my $in_obj = Bio::SeqIO->new(-file => $outfile, '-format' =>'fasta');
my $total=0;
# grab sequence object
while (my $seq = $in_obj->next_seq() ) {
	my $seq_obj = $in_obj;
	my $sequenceid = $seq->id;
	my $species_name = $seq->desc;
	my $fullheader = $sequenceid." ".$species_name;
	my $sequence = $seq->seq;
	my @header = split(/_/, $fullheader);
	my $cluster = $header[0];
	my $seqgi = $header[1];
	$seqgi =~ s/gi//;
	my $seqti = $header[2];
	$seqti =~ s/ti//;
	my $seqsp = $header[3];
	$seqsp = cleansp($seqsp);
	print "Writing phytab for $seqsp\n";
	print PHYTAB $seqsp."\t".$cluster."\t".$seqgi."\t".$sequence."\n";
}
close(PHYTAB);

#remove duplicate lines (trees)
sub uniq {
    my %seen = ();
    my @r = ();
    foreach my $a (@_) {
        unless ($seen{$a}) {
            push @r, $a;
            $seen{$a} = 1;
        }
    }
    return @r;
}
sub cleansp
{
	my $seqsp = shift;
	$seqsp =~ s/ /_/g;
	$seqsp =~ s/\.//g;
	$seqsp =~ s/\'//g;
	$seqsp =~ s/\-//g;
	return($seqsp);
}
sub getfastafromphylota
{
	my $ci=shift;
	my $ti=shift;

	print "Writing: CI:$ci TI:$ti\n";

	my $url = 'http://phylota.net/cgi-bin/sql_getcluster_fasta.cgi?format=all&db=184&ti='.$ti.'&cl='.$ci.'&ntype=1';
	my $content = get $url;
	die "Couldn't get $url" unless defined $content;
	$content =~ s/\<html\>\<pre\>//;
	$content =~ s/\<\/html\>//;
	$content =~ s/\<\/pre\>//;
	return($content);
}
