#! /usr/bin/perl -w

use strict;
use warnings;

#Written by Todd H. Oakley UCSB
#



#Obtain Arguments
my $infile1=shift(@ARGV);		#0 input file
my $minsp_pergf=shift(@ARGV);		#1 minimum number of genes to keep species
my $min_g_len =shift(@ARGV);		#3 minimum gene family length
my $outfile=shift(@ARGV);		#6 data outfile
my $logfile=shift(@ARGV);		#8 log file


my $ignoregaps = 1;

#Open 2 output files
open (OUT, ">$outfile") or die "Cannot create $outfile \n";
open (LOG, ">$logfile") or die "Cannot create $logfile\n";

my %HoAdatatable;
my %HoGF;
my %lengthof;
my %ModelFor;
my $line = "";

my @specieslist;
my @genelist;
my @nsp;
my $nsl=1;
my @uniquemodels;
my @nullmodels;
my @retainedspecies;


#unused
my $modelsfile;
my $mingf_persp;
my $speciesfile;

#INPUT ALL DATA INTO HASH
open(INFILE1, $infile1);
foreach $line(<INFILE1>) {
	chop($line);
	my $getline = $line;
	my @column = split(/\t/, $getline);
	my $species = $column[0];
	my $genefamily = $column[1];
	my $genename = $column[2];
	my $sequence = $column[3];

	if (exists $HoAdatatable{$species}{$genefamily}) {
		print OUT "ERROR: $species $genefamily is duplicated\n";
		die;
	}else{
		my $seqlen;
		if($ignoregaps)	{
			my $ungapped = $sequence;
			if(defined($ungapped)){
				$ungapped =~ s/\-//g;
				$seqlen = length($ungapped);
			}else{
				$seqlen=0;
			}
		}else{
			$seqlen = length($sequence);
		}
		if($seqlen > $min_g_len){	#sequence is long enough so keep it
			$HoGF{$genefamily}{$species}=$genename."\t".$sequence;
		}else{
			print LOG "$species \t $genefamily \t $genename \t TOO SHORT\n";
		}			
		push(@specieslist,$species);
	}
}
@specieslist = uniq(@specieslist); 	#uniq subroutine defined at end

foreach my $gfkey (sort keys %HoGF) {
	my $numsp = 0;
	for(my $j=0; $j<@specieslist; $j++){
                if(exists $HoGF{$gfkey}{$specieslist[$j]}){
			$numsp ++;
		}
	}
	if($numsp >  $minsp_pergf) { 
		for(my $j=0; $j<@specieslist; $j++){
	                if(exists $HoGF{$gfkey}{$specieslist[$j]}){
				print OUT $specieslist[$j]."\t".$gfkey."\t".$HoGF{$gfkey}{$specieslist[$j]};
				print OUT "\n";
			}
		}
	}else{
		print LOG "Deleting Gene Family $gfkey. Only Present from $numsp species\n";
	}
}
exit;

sub uniq {
    return keys %{{ map { $_ => 1 } @_ }};
}
