#!/usr/bin/perl

my $input = $ARGV[0];
my $excludeList = $ARGV[1];

my @species = split(/,/, $excludeList);

open(INPUT, $input);
	open(OUTPUT, '>'."output.txt");
		while(my $currLine = <INPUT>) {
			my @currentLine = split(/\t/, $currLine);
			my $flag = 0;
			for(my $i = 0; $i < @species; $i++) {
				if($species[$i] eq $currentLine[0]) {
					$flag = 1;
				}
			}
			if($flag == 0) {
				print OUTPUT $currLine;
			}			
		}
	close(OUTPUT);
close(INPUT);