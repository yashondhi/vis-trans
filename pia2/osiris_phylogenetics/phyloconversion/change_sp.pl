#!/usr/bin/perl -w
use strict;

my $infile = $ARGV[0];
my $changefile = $ARGV[1];
my $outfile = $ARGV[2];

open(IN, "$infile") or exit;
open(CHANGE, "$changefile");
open(OUT, ">$outfile") or exit;


my %speciesFor;	#Hash to associate code with species name

while (<CHANGE>) 
{
	chomp;
	my $currentinput = "$_";
	if($currentinput =~m /\t/){     #must have a tab otherwise wrong file format
		if($currentinput =~m /\t\t/){
			print OUT "ERROR: file contains 2 tabs in a row.  Check phytab format.\n";
                        die("ERROR: file contains 2 tabs in a row. Check it is in phytab format");
                }else{
                        my @changepair = split(/\t/, $currentinput);
                        my $codename=$changepair[0];
                        my $sp_name = $changepair[1];
                        if (exists $speciesFor{$codename}) {
                        	print OUT "ERROR: Species name specification for $codename is duplicated\n";
                                die("ERROR: Species name specifiation for for $codename is duplicated\n");
                        }else{
                        	$speciesFor{$codename}=$sp_name;
			}
                }
        }else{
		die "ERROR: Species conversion table  must be genefamily\tmodel and contain no blank lines\n";
        }
}
while (<IN>) {
	chomp;
	my $currentinput = "$_";
	if($currentinput =~m /\t/){     #must have a tab otherwise wrong file format
		if($currentinput =~m /\t\t/){
			print OUT "ERROR: file contains 2 tabs in a row.  Check phytab format.\n";
                        die;
                }else{
                        my @changepair = split(/\t/, $currentinput);
                        my $sp_name=$changepair[0];

                        if (exists $speciesFor{$sp_name}) {
				$currentinput =~s /$sp_name/$speciesFor{$sp_name}/ ;
                        	print OUT $currentinput."\n";
                        }else{
				print OUT $currentinput."\n";
			}
                }
        }else{
		die "ERROR: Input a PHYTAB file in Tabular format\n";
        }
}
close(IN);
close(OUT);
close(CHANGE);

sub change
{

	my $changetext = shift;
	$changetext =~ s/ /_/g;
	return $changetext;
}
