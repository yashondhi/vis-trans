#!/usr/bin/perl

my $file1 = $ARGV[0];
my $file2 = $ARGV[1];

my $run = qx/phylomatic -f $file1 -t $file2 > output.txt 2> errors.txt /;

print $run;
