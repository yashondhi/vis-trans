#! /usr/bin/perl -w

use strict;
use warnings;

#Written by Todd H. Oakley UCSB


#Obtain Arguments
my $table=shift(@ARGV);			#0 input file
my $infile1=shift(@ARGV);		#0 input file
my $hmmfile=shift(@ARGV);		#1 hmm definition file
my $outfile=shift(@ARGV);		#1 outfile of hmmsearch
my $outfile2=shift(@ARGV);		#2 manipulated table outfile

system "hmmsearch --$table $outfile $hmmfile $infile1";

my $finalcol;
if($table eq "tblout"){
	$finalcol=18;
}elsif($table eq "domtblout"){
	$finalcol=22;
}	

#READ in table output from hmmsearch, then recast as tab delimited
open(HMMSEARCHFILE, $outfile);
open(OUT, ">$outfile2");
foreach my $line(<HMMSEARCHFILE>) {
	chop($line);
	my $getline = $line;
	if($getline =~ m/^# +\-/ ){
		next;
	}
	$getline =~ s/^\# /#/;
	$getline =~ s/target name/target_name/;
	$getline =~ s/query name/query_name/;
	$getline =~ s/description of target/description_of_target/;
	my @column = split(/ +/, $getline);
	print OUT join("\t", @column[0..$finalcol-1]);
	print OUT "\t";
	if(@column > $finalcol){
		print OUT join(" ", @column[$finalcol..$#column]);
	}else{
		print OUT $column[$finalcol];
	}
	print OUT "\n"; 
}
close(HMMSEARCHFILE);
close(OUT);
