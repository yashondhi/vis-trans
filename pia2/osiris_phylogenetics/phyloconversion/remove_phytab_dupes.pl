#!/usr/bin/perl
use strict;



my $infile=$ARGV[0];
my $keeplongest=$ARGV[1];
my $ignoregaps=$ARGV[2];
my $uniout=$ARGV[3];
my $dupout=$ARGV[4];

open IN, $infile or die "Cannot open $infile\n";

my %UniquesHash;
my @DupeArray;

while(<IN>){
	my $row = $_;
	chomp($row);
	my @column = split(/\t/, $row);
	my $species = $column[0];
	my $partition = $column[1];
	my $id = $column[2];
	my $sequence = $column[3];

	if(exists $UniquesHash{$species}{$partition}){
		my @dupeseq = split(/\t/, $UniquesHash{$species}{$partition});
		my ($savlen,$curlen);
		if($ignoregaps==1){
			my $nogapsav = $dupeseq[1];
			my $nogapcur = $sequence;
			$nogapsav =~ s/\-//g;
			$nogapcur =~ s/\-//g;
			$savlen = length($nogapsav);
			$curlen = length($nogapcur);
		}else{
			$savlen = length($dupeseq[1]);
			$curlen = length($sequence);
		}
		if($curlen > $savlen && $keeplongest==1) { 		#current is longer so keep that one
			my $oldline = $species."\t".$partition."\t".$UniquesHash{$species}{$partition}."\n";
			$UniquesHash{$species}{$partition} = "$id\t$sequence";
			push(@DupeArray, $oldline);
		}else{ 
			push(@DupeArray, "$species\t$partition\t$id\t$sequence\n");
		}
	}else{
		$UniquesHash{$species}{$partition} = "$id\t$sequence";
	}
}

open OUT, ">".$uniout or die "Cannot open $uniout\n";
open DUPES, ">".$dupout or die "Cannot open $dupout\n";

print DUPES @DupeArray;
for my $spname ( keys %UniquesHash ) {
    for my $partname ( keys %{ $UniquesHash{$spname} } ) {
         print OUT "$spname\t$partname\t$UniquesHash{$spname}{$partname}\n";
    }
}
