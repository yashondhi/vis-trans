#!/usr/bin/perl

my $file = $ARGV[0];
my $q_mark = $ARGV[1];
my $hyphen = $ARGV[2];
my $N = $ARGV[3];
my $usePartFile = $ARGV[4];
my $partFile = $ARGV[5];

my $out = "out.phylipnon";	# output file

open(FILE, $file);
	my @speciesNames;
	my @sequenceLines;
	
	my @currentLineContent;
	
	my $i = 0;
	while($currentLine = <FILE>) {
		chomp($currentLine);
		@currentLineContent = split(/\t/, $currentLine);
		$speciesNames[$i] = $currentLineContent[0];
		$sequenceLines[$i] = $currentLineContent[1];
		$i++;
	}
	
	my $dataInfo = $speciesNames[1];	# gets num of species and sequence length
	my @numbers = split(/ /, $dataInfo);

	my $numberOfSpecies = $numbers[0];
	my $sequenceLength = $numbers[1];
	
close(FILE);

open(OUT, '>'.$out);
	my @columnData;		# this will have $sequenceLength elements
	for($j = 0; $j < $numberOfSpecies+2; $j++) {
		for($k = 0; $k < $sequenceLength; $k++) {
			$currChar = substr($sequenceLines[$j], $k, 1);
			$columnData[$k] = $columnData[$k].$currChar;
		}
	}
	
	# mark locations that will be removed
	my @flagMap;
	for($i = 0; $i < $sequenceLength; $i++) {
		$flagMap[$i] = 0;		
	}
	my $index = 0;
	foreach $el(@columnData) {
		my $tot = 0;
		my $q_mark_occur = 0;
		my $hyphen_occur = 0;
		my $N_occur = 0;
		
		if($q_mark eq "true") {
			$q_mark_occur = ($el =~ tr/?//);
		}
		if($hyphen eq "true") {
			$hyphen_occur = ($el =~ tr/-//);	
		}
		if($N eq "true") {
			$N_occur = ($el =~ tr/N//);
		}

		$tot = $q_mark_occur + $hyphen_occur + $N_occur;
		if($tot == $numberOfSpecies) {
			$flagMap[$index] = 1;
		}
		$index++;
	}
	
	my $newSequenceLength = $sequenceLength;
	foreach $el(@flagMap) {
		if($el == 1) {
			$newSequenceLength--;
		}
	}

	print OUT $speciesNames[0]."\n";
	print OUT $numberOfSpecies." ".$newSequenceLength."\n";
	for($i = 2; $i < $numberOfSpecies+3; $i++) {
		print OUT $speciesNames[$i]."\t";
		for($j = 0; $j < $sequenceLength; $j++) {
			if($flagMap[$j] == 0) {
				my $character = substr($sequenceLines[$i], $j, 1);
				print OUT $character;
			}
		}
		print OUT "\n"; 
	}	

close(OUT);

my $partOut = "partOut.txt";

if($usePartFile eq "true") {
	# update the partition file
	open(PART, $partFile);
		my @data;
		my @ranges;
		my @names;
		$i = 0;
		while($currentLine = <PART>) {
			@data = split(/=/, $currentLine);
			$names[$i] = $data[0];
			$ranges[$i] = $data[1];
			$i++;
		}
	close(PART);
	
	my $firstFlag = 1;
	open(PARTOUT, '>'.$partOut);
		$j = 0;
		my $newLower;
		foreach $el(@ranges) {
			print PARTOUT $names[$j]." = ";
			@lowerUpper = split(/-/, $el);
			if($firstFlag == 1) {
				$newLower = $lowerUpper[0];
				$firstFlag = 0;
			}
			my $currUpper = $lowerUpper[1];	
			my $newUpper = $currUpper;

			

			for($i = $currLower; $i < $currUpper; $i++) {
				if($flagMap[$i] == 1) {
					$newUpper--;
				}
			}

			print PARTOUT $newLower." - ".$newUpper."\n";
			$newLower = $newUpper + 1;
			$j++;
		}
	close(PARTOUT);
}
