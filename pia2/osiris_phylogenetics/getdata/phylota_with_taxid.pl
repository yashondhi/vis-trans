#!/usr/bin/env perl
use strict;
use LWP::Simple;
use Bio::SeqIO;

my $ti = $ARGV[0];
my $outfile = $ARGV[1];
my $phytabfile = $ARGV[2];

open(OUT, ">$outfile") or exit;

my $content = getclustersfromphylota($ti);
my @weblines = split(/\<\/tr\>/, $content);
my @ci;

#Parse html from phylota browser to retain just each ci
foreach(@weblines){
	if($_ =~ m/getcluster\.cgi/){
		chomp;
		$_ =~ s/\&ntype\=1\&db\=184\".+// ;
		$_ =~ s/(.*?)getcluster\.cgi.+cl\=// ;
		$_ =~ s/\<\/font\>\<\/td\>// ;
		chomp;
		$_ =~ s/^\n// ;
		push(@ci, $_);
	}
}

#get fasta files for trees
for(my $i=0;$i < @ci; $i++){
	my $ci = $ci[$i];
	my $addstring = 'ti'.$ti.'ci'.$ci.'_';
	my $fastafile = getfastafromphylota($ci,$ti);
	#Add TI_CI_ to each fastaheader
	$fastafile =~ s/\>/\>$addstring/g;
	print OUT $fastafile;
}
close(OUT);

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
	print PHYTAB $seqsp."\t".$cluster."\t".$seqgi."\t".$sequence."\n";
}
close(PHYTAB);






#**************************************************************
#sub routines

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

	#print "Writing: CI:$ci TI:$ti\n";

	my $url = 'http://phylota.net/cgi-bin/sql_getcluster_fasta.cgi?format=all&db=184&ti='.$ti.'&cl='.$ci.'&ntype=1';
	my $content = get $url;
	die "Couldn't get $url" unless defined $content;
	$content =~ s/\<html\>\<pre\>//;
	$content =~ s/\<\/html\>//;
	$content =~ s/\<\/pre\>//;
	return($content);
}
sub getclustersfromphylota
{
	my $ti=shift;

	#print "Writing: CI:$ci TI:$ti\n";

	my $url = 'http://phylota.net/cgi-bin/sql_getclusterset.cgi?ti='.$ti.'&ntype=1&piflag=1&dflag=0&db=184';
	my $content = get $url;
	die "Couldn't get $url" unless defined $content;
	$content =~ s/\<html\>\<pre\>//;
	$content =~ s/\<\/html\>//;
	$content =~ s/\<\/pre\>//;
	return($content);
}
