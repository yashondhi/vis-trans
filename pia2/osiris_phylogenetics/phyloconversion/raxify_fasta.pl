#!/usr/bin/perl -w
use strict;

my $infile = $ARGV[0];
my $outfile = $ARGV[1];

open(IN, "$infile") or exit;
open(OUT, ">$outfile") or exit;

while (<IN>) {
	my $line = $_;
	if($line =~ m/\>/ ){
		$line=raxify($line);
	}
	if($line =~ m/\n/){
		print OUT $line;
	}else{
		print OUT $line."\n";
	}
}
close(IN);
close(OUT);

sub raxify
{
	my $raxline = shift;
	$raxline = substr($raxline,0,51);
	$raxline =~ s/\./_/g;
	$raxline =~ s/\|/_/g;
	$raxline =~ s/ /_/g;
	return $raxline;
}
