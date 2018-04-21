#!/usr/bin/perl -w
use strict;

my $datafile = $ARGV[0];
my $keepfile = $ARGV[1];
my $delfile = $ARGV[2];
my $subsp = $ARGV[3];
my $var = $ARGV[4];

open (FILE,"<$datafile") or die "Cannot open file input file\n";
open (KFILE,">$keepfile") or die "Cannot open file $keepfile\n";
open (DFILE,">$delfile") or die "Cannot open file delfile\n";

my $keep = 1;

while (<FILE>)
{
	if($_ =~ m/_\d/){
		$keep=0;
	}else{
		if($subsp==1){
			if($_ =~ m/subsp/){
				$keep=0;
			}
		}
		if($var==1){
			if($_ =~ m/_var_/){
				$keep=0;
			}
		}
	}

	if($keep == 0){
		print DFILE $_;
	}else{
		print KFILE $_;
	}
	$keep=1; #reset variable. Default is keep
}
close FILE;
close KFILE;
close DFILE;

