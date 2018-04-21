#!/usr/bin/env perl
#
# Modified for Osiris by THO.
# 
#
# GenBankStrip.pl v2.0
# Last modified July 25, 2005 14:29
# (c) Olaf R.P. Bininda-Emonds
#
# Input:
#   1) A GenBank output file
#   2) A file containing gene names, each on a separate line
#
# Generates:
#   1) A Se-Al formatted (and optionally a nexus-formatted) datafile containing stripped sequences
#	2) A summary file giving status of each entry in GenBank file
#
# Usage: GenBankStrip.pl -f<filename> [-g<filename>] [-k<number>] [-l<number>] [-o<n|s>] [-s] [-t] [-h] [-v]
#   options: -f<filename> = file containing sequences in GenBank format
#            -g<filename> = file containing specific genes to be stripped
#            -i<genename> = strip a single gene
#            -k<number> = number of (longest) sequences to retain per species for a given gene (default = all)
#            -l<number> = minimum length required for all non-tRNA genes (default = none)
#            -o<n|s> = provide output in nexus (n) phytab (p) and/or Se-Al (s) format in addition to fasta format
#            -s = only process sequences with valid species names (i.e., no species with sp. or cf. in names)
#            -t<g<number>|s<number>> = minimum number of sequences (g; default = 0) or species (s; default = 0) a gene must have to be retained; both can be specified simultaneously
#            -h = print this message and quit
#            -v = verbose output
#	     -a = Accession number appears first in stripped Fasta File -- Added by THO
use strict;

# Set default values
	# Genes and sequences
		my %synonym;
		my %complement = ('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A',
						 'M' => 'K', 'R' => 'Y', 'W' => 'S', 'S' => 'W', 'Y' => 'R', 'K' => 'M',
						 'B' => 'V', 'D' => 'H', 'H' => 'D', 'V' => 'B', 'N' => 'N', '-' => '-', '?' => '?');
		my %maxLength;
	
	# I/O
		my $gbFile;
			my ($countFile, $rejectFile, $stripFile);
		my $geneFile;
			my %userGene;
		my $singleGene;
		my $keepGene = 0;
		my $minLength;
			my $minLengthTRNA;
		my $speciesBlock = 0;
		my $seqThreshold = 0;
		my $speciesThreshold = 0;
		my $sealPrint = 0;
		my $nexusPrint = 0;
		my $phytabPrint = 0;
		
	# Global stats
		my (@globalGeneList, %globalGenePresent, %geneStatus, @rejectList, $stripCount);
		my (@speciesList, %speciesPresent, %speciesGenePresent, %quickSpeciesCount, %speciesCount, %quickSequenceCount, %sequenceCount);
				$sequenceCount{"all"} = 0;

	# Miscellaneous
		my $debug = 0;
		my $verbose = 0;
		my $accessionFirst = 0; #Added by THO as option to have accession # first in FASTA output

for (my $i = 0; $i <= $#ARGV; $i++)
	{
#Original stripped out directory structure, required by Galaxy
#	if ($ARGV[$i] =~ /^-f([\w+\.?\w+?]*)/)
	if ($ARGV[$i] =~ /^-f(.+)/)
		{
		$gbFile = $1;
		
#		(my $baseFile = $gbFile) =~ s/.\w+$//;
#			$countFile .= $baseFile . "_geneCount.txt";
#			$rejectFile .= $baseFile . "_gbs.rejectlist.txt";
#			$stripFile .= $baseFile . "_gbs.striplist.txt";
#		}
#for Galaxy, easier to have same file name each time it'r run
		(my $baseFile = $gbFile) =~ s/.\w+$//;
			$countFile .= "geneCount.txt";
			$rejectFile .= "rejectlist.txt";
			$stripFile .= "striplist.txt";
		}
	elsif ($ARGV[$i] =~ /^-g([\w+\.?\w+?]*)/)
		{
		$geneFile = $1;

#This is a hack to write the name of a single gene into a new genefile based on the Galaxy call

open (GENE,">$geneFile") or die "Cannot open file $geneFile containing gene names.\n";
print GENE "$geneFile\n";
close GENE;
		}
	elsif ($ARGV[$i] =~ /^-i([\w+\.?\w+?]*)/)
		{
		$singleGene = $1;
		}
	elsif ($ARGV[$i] =~ /^-k(\d+)/)
		{
		$keepGene = $1;
		}
	elsif ($ARGV[$i] =~ /^-l(\d+)/)
		{
		$minLength = $1;
		$minLengthTRNA = 50;
		}
	elsif ($ARGV[$i] =~ /^-on/)
		{
		$nexusPrint = 1;
		}
	elsif ($ARGV[$i] =~ /^-op/)
		{
		$phytabPrint = 1;
		}
	elsif ($ARGV[$i] =~ /^-os/)
		{
		$sealPrint = 1;
		}
	elsif ($ARGV[$i] =~ /^-s/)
		{
		$speciesBlock = 1;
		}
	elsif ($ARGV[$i] =~ /^-tg(\d+)/)
		{
		$seqThreshold = $1;
		}
	elsif ($ARGV[$i] =~ /^-ts(\d+)/)
		{
		$speciesThreshold = $1;
		}
	elsif ($ARGV[$i] =~ /^-v/)
		{
		$verbose = 1;
		}
	elsif ($ARGV[$i] =~ /^-a/)
		{
		$accessionFirst = 1;
		}
	elsif ($ARGV[$i] =~ /^-x/)
		{
		$debug = 1;
		$verbose = 1;
		}
	elsif ($ARGV[$i] =~ /^-h/)
		{
		print "Usage: GenBankStrip.pl -f<filename> [-g<filename>] [-k<number>] [-l<number>] [-o<n|s>] [-s] [-t] [-h] [-v]\n";
		print "Options: -f<filename> = file containing sequences in GenBank format\n";
		print "         -g<filename> = file containing specific genes to be stripped\n";
		print "         -k<number> = number of (longest) sequences to retain per species for a given gene (default = all)\n";
		print "         -l<number> = minimum length required for all non-tRNA genes (default = none)\n";
		print "         -o<n|s> = provide output in nexus (n) phytab (p) and/or Se-Al (s) format in addition to fasta format\n";
		print "         -s = only process sequences for valid species names (i.e., no species with sp. or cf. in names)\n";
		print "         -t<g<number>|s<number>> = minimum number of sequences (g; default = 0) or species (s; default = 0)\n";
		print "                                   a gene must have to be retained; both can be specified simultaneously\n";
		print "         -h = print this message and quit\n";
		print "         -v = verbose output\n";
		exit(0);
		}
	else
		{
		print "Don't understand argument: $ARGV[$i]\n";
		print "Usage: GenBankStrip.pl -f<filename> [-g<filename>] [-k<number>] [-l<number>] [-o<n|s>] [-s] [-t] [-h] [-v]\n";
		exit(1); 
		}
	}

die "ERROR: Must supply name of GenBank output file using flag -f.\n" if (not $gbFile);

# Load in hardwired gene synonyms
	geneSynonyms();

# Get list of target genes (if desired)
	if ($geneFile)
		{
		my $userGeneCount = 0;
		my %userGenePresent;
		setLineBreak($geneFile);
		open (GENE,"<$geneFile") or die "Cannot open file $geneFile containing gene names.\n";
		print "Gene(s) to be stripped:\n";
		while (<GENE>)
			{
			chomp;
			next unless ($_);
			my $gene = $_;
			$gene = $synonym{$gene} if (defined $synonym{$gene});
			$userGene{$gene} = 1;
			unless ($userGenePresent{$gene})
				{
				$userGeneCount++;
				$userGenePresent{$gene}++;
				print "\t$gene\n";
				}
			}
		close GENE;

#		die "ERROR: No genes read in from file $geneFile\n";
		}
#THO added -h command to easilly pass a single gene for Galaxy
	if($singleGene)
		{
		my $userGeneCount = 1;
		my %userGenePresent;
		print "Gene(s) to be stripped:\n";
		my $gene = $singleGene;
		$gene = $synonym{$gene} if (defined $synonym{$gene});
		print "\t$gene\n";
		$geneFile = $singleGene;
		}

# Print parameter summary
	print "The following parameters have been set by the user:\n";
	print "\tFile containing GenBank sequence data: $gbFile\n";
	print "\tFile containing target genes to be stripped: $geneFile\n" if ($geneFile);
	print "\tUser-defined constraints(s):\n";
		print "\t\tNumber of sequences: $seqThreshold\n" if ($seqThreshold);
		print "\t\tNumber of species: $speciesThreshold\n" if ($speciesThreshold);
		print "\t\tMinimum sequence length: global - $minLength bp; tRNAs - $minLengthTRNA bp\n" if (defined $minLength);
		print "\t\tOnly using species with valid names\n" if ($speciesBlock);
		print "\t\tNone\n" if (not $seqThreshold and not $speciesThreshold and not defined $minLength and not $speciesBlock);
	print "\tNumber of sequences to keep per species for each gene: ";
		if ($keepGene)
			{
			print "$keepGene\n";
			}
		else
			{
			print "all\n";
			}
	print "\tOutput format(s): fasta";
		print ", Se-Al" if ($sealPrint);
		print ", nexus" if ($nexusPrint);
		print ", phytab" if ($phytabPrint);
		print "\n";

# Do quick gene count if thresholds indicated; takes longer but 1) saves many disk operations and 2) less memory intensive
	geneCount($gbFile) if (not defined $geneFile and ($seqThreshold or $speciesThreshold));	# Don't bother if user gene list given; will usually be small enough so that benefits won't come into play
	
# Read in GenBank file and strip genes
	my $stripZero = time;
	print "\nProcessing GenBank file $gbFile ...\n";
	setLineBreak($gbFile);
	open (DATA, "<$gbFile") or die "Cannot open GenBank output file $gbFile\n";
		my @allAccNum;
		my ($accNum, %accRead, $duplEntry, $organism, %species, $geneName);
			my $speciesFlag = 0;
		my (%genePresent, @geneList);
		my (@startList, @stopList, $complementFlag, $typeStatus, %pseudoStatus, %seqType, $joinLine, @joinSegments, %geneStart, %geneStop, %compStatus, $fullSeq);
			my $nameFlag = 0;
			my $joinFlag = 0;
			my $pseudoFlag = 0;
			my $readSeqs = 0;

		while (<DATA>)
			{
			chomp;
			my $readLine = $_;
			next if ($readLine =~ /^\s*LOCUS/ or $readLine =~ /^\s*DEFINITION/ or $readLine =~ /^\s*VERSION/ or $readLine =~ /^\s*KEYWORDS/ or $readLine =~ /^\s*SOURCE/ or $readLine =~ /^\s*ORGANISM/ or $readLine =~ /^\s*REFERENCE/ or $readLine =~ /^\s*AUTHORS/ or $readLine =~ /^\s*TITLE/ or $readLine =~ /^\s*JOURNAL/ or $readLine =~ /^\s*MEDLINE/ or $readLine =~ /^\s*PUBMED/ or $readLine =~ /^\s*FEATURES/ or $readLine =~ /^\s*COMMENT/ or $readLine =~ /^\s*BASE COUNT/ or $readLine =~ /^\s*source/ or $readLine =~ /^\s*\/codon/ or $readLine =~ /^\s*\/transl/ or $readLine =~ /^\s*\/db_/ or $readLine =~ /^\s*CONTIG/);
		
			# Get accession number
				if ($readLine =~ /^\s*ACCESSION\s+(.+)/)
					{
					$accNum = $1;
					print "$readLine\n" if ($debug);
					# Clear variables
						undef @geneList;
						undef %genePresent;
						undef $fullSeq;
						undef @startList;
						undef @stopList;
						undef %geneStart;
						undef %geneStop;
						undef %pseudoStatus;
						undef $geneName;
						$speciesFlag = $nameFlag = $joinFlag = $pseudoFlag = $readSeqs = $duplEntry = 0;

					if (not $accRead{$accNum})	# Check for duplicate entries
						{
						$accRead{$accNum} = 1;
						push @allAccNum, $accNum;
						if (scalar(@allAccNum) == int(scalar(@allAccNum)/10000)*10000 and $verbose)
							{
							print "\tSequences read in: ".scalar(@allAccNum)."\n";
							}
						print "\tAccession number: $accNum\n" if ($debug);
						}
					else
						{
print "*****NOTE --  DUPLICATE ENTRY Accession $accNum\n";
						$duplEntry = 1;
						}
					}

			# Get organism name
				if ($readLine =~ /^\s*\/organism=\"(.+)\"/)
					{
					$organism = $1;
						$organism =~ s/\s/_/g;
					print "$readLine\n" if ($debug);
					$species{$accNum} = $organism;
					$speciesFlag = 1 if ($organism =~ /sp\./ or $organism =~ /cf\./ or $organism =~ /_X_/i);
					if ($debug)
						{
						print "\t\tOrganism: $organism";
						print " (blocked)" if ($speciesFlag and $speciesBlock);
						print "\n";
						}
					}
				next if ($speciesFlag and $speciesBlock);	# Entry pertains to invalid species name; skip parsing rest of entry

			
			# Get gene boundaries; process previous set of CDs
				if ($readLine =~ /\<?(\d+)\<?\.\.\>?(\d+)\>?/ or $joinFlag == 1 or $readLine =~ /^\s*ORIGIN/)	# ORIGIN will process last set of CDs
					{
					next if ($readLine =~ /^\s+\/\w+/);	# Prevents spurious matches with lines beginning with "/feature"
					$readSeqs = 1 if ($readLine =~ /^\s*ORIGIN/);	# Indicates that remaining lines will contain sequence information
					
					# Process previous gene; need to do here to account for a posteriori declarations of pseudogene status
						if ($geneName and $nameFlag == 2 and @startList and @stopList)	# Process complete gene
							{
							print "\t\t\t\tParsed name: $geneName (type: $typeStatus)\n" if ($debug);
							
							# Clean up gene name and misleading punctuation
								$geneName = geneClean($geneName);
								$pseudoStatus{$geneName} = 1 if ($pseudoFlag);
								
							if (defined @ { $geneStart{$geneName} } and ((defined $pseudoStatus{$geneName} and $pseudoStatus{$geneName} == 1) or $typeStatus =~ /intron/i or $typeStatus =~ /UTR/i))	# Gene has previously stored CDs that might not have been recognized as a pseudogene or non-coding region
								{
								print "\t\t\t\t\tSubsequent occurrence of now recognized pseudogene or non-coding region; comparing new and stored CDs\n" if ($debug);
								for (my $i = 0; $i < scalar(@startList); $i++)	# Check each occurrence in new CDs for matches in stored ones
									{
									my $newStart = $startList[$i];
									my $newStop = $stopList[$i];
									print "\t\t\t\t\t\tChecking new CDs $newStart to $newStop\n" if ($debug);
									for (my $j = 0; $j < scalar(@ { $geneStart{$geneName} }); $j++)
										{
										if ($newStart == $geneStart{$geneName}[$j] and $newStop == $geneStop{$geneName}[$j])
											{
											print "\t\t\t\t\t\t\tMatch with stored CDs (no. $j); deleted\n" if ($debug);
											splice(@ { $geneStart{$geneName} }, $j, 1);
											splice(@ { $geneStop{$geneName} }, $j, 1);
											}
										}
									}
								if ($debug)
									{
									print "\n\t\t\t\t\tCurrent gene boundaries after pseudogene / non-coding check (type = $seqType{$geneName}):\n";
									if (scalar(@ { $geneStart{$geneName} }) < 1)
										{
										print "\t\t\t\t\t\tnone\n";
										}
									else
										{
										for (my $i = 0; $i < scalar(@ { $geneStart{$geneName} }); $i++)
											{
											print "\t\t\t\t\t\t$geneStart{$geneName}[$i]\t$geneStop{$geneName}[$i]\n";
											}
										}
									}
								}
														
							# Only process coding regions of user-desired genes (if applicable), genes with sensible CDs, non-blocked genes, and genes that are not obvious singletons or pseudogenes
								unless (($geneFile and not defined $userGene{$geneName}) or
									(defined $geneStatus{$geneName} and $geneStatus{$geneName} eq "rejected") or
									(defined $pseudoStatus{$geneName} and $pseudoStatus{$geneName} == 1) or
									$geneName =~ /hypothetical/i or $geneName =~ /open reading frame/i or $geneName =~ /similar/i or $geneName =~ /homolog/i or $geneName =~ /putative/i or $geneName =~ /unknown/i or $geneName =~ /unnamed/i or $geneName =~ /\d+rik/i or $geneName =~ /possible/i or $geneName =~ /pseudo/i
									or $typeStatus =~ /UTR/i or $typeStatus =~ /intron/i or $typeStatus =~ /misc_feature/i)
									{
									if (not defined $genePresent{$geneName})	# Process first occurrence of gene in entry
										{
										if ($debug)
											{
											print "\t\t\t\t\tFirst occurrence of $geneName\n" if ($debug);
											for (my $i = 0; $i < scalar(@startList); $i++)
												{
												print "\t\t\t\t\t\t$startList[$i]\t$stopList[$i]\n";
												}
											}
										$genePresent{$geneName} = 1;
										push @geneList, $geneName;
											push @ { $geneStart{$geneName} }, $_ foreach (@startList);
											push @ { $geneStop{$geneName} }, $_ foreach (@stopList);
			
										$seqType{$geneName} = $typeStatus;
										$compStatus{$geneName} = 0;	# Note whether gene is complemented or not
											$compStatus{$geneName} = 1 if ($complementFlag);
										}
									else	# Attempt to add secondary occurrences
										{
										my $storedSegments = scalar(@ { $geneStop{$geneName} }) - 1;
										my $lastStop = $geneStop{$geneName}[$storedSegments];
											$lastStop = 0 if (not defined $lastStop);
										my $newStart = $startList[0];
			
										if ($debug)
											{
											print "\t\t\t\t\tSecondary occurrence of $geneName with boundaries\n" if ($debug);
											for (my $i = 0; $i < scalar(@startList); $i++)
												{
												print "\t\t\t\t\t\t$startList[$i]\t$stopList[$i]\n";
												}
											}
										if ($seqType{$geneName} eq "gene" and $typeStatus ne "gene")	# New information probably more precise and better accounts for structure
											{
											print "\n\t\t\t\t\t\tNew segment more precisely defined by type; replaced\n" if ($debug);
											undef @ { $geneStart{$geneName} };
											undef @ { $geneStop{$geneName} };
											push @ { $geneStart{$geneName} }, $_ foreach (@startList);
											push @ { $geneStop{$geneName} }, $_ foreach (@stopList);
											$seqType{$geneName} = $typeStatus;
											}
										
										elsif ($newStart > $lastStop)	# New segment occurs distal to last stored segment (could also be contiguous); append to boundaries
											{
											print "\n\t\t\t\t\t\tContiguous with or occurs after last stored segment; appended\n" if ($debug);
											push @ { $geneStart{$geneName} }, $_ foreach (@startList);
											push @ { $geneStop{$geneName} }, $_ foreach (@stopList);
											$seqType{$geneName} = "composite";
											}
										elsif (scalar(@ { $geneStart{$geneName} }) == 1 and scalar(@startList) > 1)	# Replace single stored segment with new segments derived from join statement; probably better accounts for intron/exon structure
											{
											print "\n\t\t\t\t\t\tNew segments subdivide single stored segment; replaced\n" if ($debug);
											undef @ { $geneStart{$geneName} };
											undef @ { $geneStop{$geneName} };
											push @ { $geneStart{$geneName} }, $_ foreach (@startList);
											push @ { $geneStop{$geneName} }, $_ foreach (@stopList);
											$seqType{$geneName} = "composite";
											}
										else
											{
											print "\n\t\t\t\t\t\tNew segment overlaps stored segments; rejected\n" if ($debug);
											}
										}
									if ($debug)
										{
										print "\n\t\t\t\t\tCurrent gene boundaries after CDS processing (type = $seqType{$geneName}):\n";
										if (scalar(@ { $geneStart{$geneName} }) < 1)
											{
											print "\t\t\t\t\t\tnone\n";
											}
										else
											{
											for (my $i = 0; $i < scalar(@ { $geneStart{$geneName} }); $i++)
												{
												print "\t\t\t\t\t\t$geneStart{$geneName}[$i]\t$geneStop{$geneName}[$i]\n";
												}
											}
										}
									}

							# Clear entries for gene
								$nameFlag = $complementFlag = 0;
								undef @startList;	# Prevents double listing of genes if CDS use both /gene and /product tags
								undef @stopList;
								undef $typeStatus;
								$pseudoStatus{$geneName} = 0;	# Clear pseudogene status for previous set of CDs
							}
						next if ($readLine =~ /^\s*ORIGIN/);	# No more CDs to worry about
					
					print "$readLine\n" if ($debug);

					# Reset gene name and information for current set of CDs
						undef $geneName;
						$nameFlag = $pseudoFlag = $complementFlag = 0;
							$complementFlag = 1 if ($readLine =~ /complement/); #DISABLED BY THO BECAUSE ONLY COMPLEMENTS AND DOESN'T REVERSE

					# Get new CD information
						if ($readLine =~ /^\s+\S+\s+join\(/i or $readLine =~ /^\s+\S+\s+order\(/)
							{
							$joinFlag = 1;
							undef @startList;
							undef @stopList;
							undef $joinLine;
							}
						if ($readLine =~ /^\s+(\S+)\s+/i)
							{
							$typeStatus = $1;
							}
						
						if ($joinFlag)
							{
							$joinLine .= $readLine;
							if ($readLine =~ /\)$/)	# Have accounted for multiline join statements; process
								{
								$joinFlag = 0;
								$complementFlag = 1 if ($joinLine =~ /complement/);
								
								# Clean up join statement
									$joinLine =~ s/complement\(//;
									$joinLine =~ s/>//g;
									$joinLine =~ s/<//g;
									$joinLine =~ s/^\s+\S+\s+join\(//;
									$joinLine =~ s/^\s+\S+\s+order\(//;
									$joinLine =~ s/\)+$//;
									$joinLine =~ s/\s+//g;
									print "\t\tJoin boundaries: $joinLine\n" if ($debug);
									if ($joinLine =~ /[A-Z]+\d+/i)	# Avoid join commands referring to separate accessions
										{
										undef $joinLine;
										print "\t\t\tERROR: Joining accessions; skipping\n" if ($debug);
										next;
										}
								
								# Derive pieces of join statement
									@joinSegments = split(/,/, $joinLine);
									foreach my $segment (@joinSegments)
										{
										my ($start, $stop) = split(/\.\./, $segment);
										# Check safeties
											next unless ($start and $stop);
											next if ($start =~ /\D/ or $stop =~ /\D/);
											next if ($start > $stop);
										push @startList, $start;
										push @stopList, $stop;									
										}
								
								undef $joinLine;
								$typeStatus = "join";
								}
							}
						else
							{
							my ($start, $stop);
							if ($readLine =~ /\<?(\d+)\<?\.\.\>?(\d+)\>?/)
								{
								$start = $1;
								$stop = $2;
								}
							next unless ($start and $stop);
							next if ($start > $stop);
							undef @startList;
								push @startList, $start;
							undef @stopList;
								push @stopList, $stop;
							print "\t\tParsed boundaries: $1\t$2\n" if ($debug);
							}
					}


			# Check for pseudogene status
				if ($readLine =~ /\s+\/pseudo/ or $readLine =~ /pseudogene/ or $readLine =~ /putative/)
					{
					print "$readLine\n" if ($debug);
					$pseudoStatus{$geneName} = 1 if (defined $geneName);
					$pseudoFlag = 1;
					print "\t\t\t\tWARNING: suspected pseudogene or of uncertain (\"putative\") status\n" if ($debug);
					}

				if ($readLine =~ /codon recognized/ and defined $geneName and $geneName =~ /trna/i)
					{
					print "$readLine\n" if ($debug);
					$pseudoStatus{$geneName} = 1;
					$pseudoFlag = 1;
					print "\t\t\t\tWARNING: tRNA for an alternative codon; ignored\n" if ($debug);
					}

			# Get gene name
				if ($readLine =~ /^\s*(\/gene)|(\/product)=\"(.+)/ and $nameFlag == 0)	# Get start of gene name
					{
					$geneName = $1 if ($readLine =~ /=\"(.+)/);
					print "$readLine\n" if ($debug);
					
					# Check whether name wraps onto a new line
						if (substr($readLine, -1) ne "\"")
							{
							$nameFlag = 1;
							next;
							}
						else
							{
							$nameFlag = 2;
							$geneName =~ s/\"$// if ($geneName =~ /\"$/);
							}
					}

				if ($nameFlag == 1)	# Get continuation / end of gene name
					{
					print "$readLine\n" if ($debug);
					$geneName .= $readLine;
					$nameFlag = 2 if ($readLine =~ /\"$/)	# Gene name is complete
					}
			
			# Read in sequence information and append
				if ($readLine =~ /^\s*\d*1\s\w+/ and @geneList and $readSeqs)
					{
					my $seqFrag = $readLine;
						$seqFrag =~ s/\s+//g;
						$seqFrag =~ s/\d+//g;
						$fullSeq .= uc($seqFrag);
					}

			# End of entry; process
				if ($readLine =~ /\/\//)
					{
					next if ($duplEntry);
					next unless (@geneList and defined $fullSeq);
					foreach my $gene (@geneList)
						{
						# Safeties; shouldn't come into play
							next if ($geneFile and not defined $userGene{$gene});
							next if (defined $geneStatus{$gene} and $geneStatus{$gene} eq "rejected");
							next if (defined $pseudoStatus{$gene} and $pseudoStatus{$gene} == 1);
						
						print "\n\t\tProcessing gene $gene\n" if ($debug);

						# Strip gene out of full sequence and format as fasta
							next unless (@ { $geneStart{$gene} } and @ { $geneStop{$gene} } );
							
							my $geneSeq;
								for (my $i = 0; $i < scalar(@ { $geneStart{$gene} }); $i++)
									{
									print "\t\t\tGetting sequence between positions $geneStart{$gene}[$i] and $geneStop{$gene}[$i]\n" if ($debug);
									next if ($geneStart{$gene}[$i] > length($fullSeq) or $geneStop{$gene}[$i] > length($fullSeq));
									my $geneSegment = substr($fullSeq,$geneStart{$gene}[$i]-1,$geneStop{$gene}[$i]-$geneStart{$gene}[$i]+1);
										print "\t\t\t\t$geneSegment\n" if ($debug);
										$geneSeq .= $geneSegment;
									}

							next unless ($geneSeq);
							$geneSeq = complement($geneSeq) if ($compStatus{$gene});
							$maxLength{$gene} = length($geneSeq) if (not defined $maxLength{$gene} or length($geneSeq) > $maxLength{$gene});

							# Check if sequence matches length threshold (if appropriate)
								if (defined $minLength)
									{
									if ($gene =~ /^trna-\w\w\w/) 
										{
										if (length($geneSeq) < $minLengthTRNA)
											{
											printf "\t\t\tRejected: length (%s bp) does not meet tRNA length threshold ($minLengthTRNA bp)\n", length($geneSeq) if ($debug);
											next;
											}
										}
									else
										{
										if (length($geneSeq) < $minLength)
											{
											printf "\t\t\tRejected: length (%s bp) does not meet global threshold ($minLength bp)\n", length($geneSeq) if ($debug);
											next;
											}
										}
									}
								if ($gene =~ /trna-\w\w\w/ and length($geneSeq) > 100)
									{
									printf "\t\t\tRejected: length of tRNA too long (%s bp); indicates parsing failure\n", length($geneSeq) if ($debug);
									next;
									}

							my $breakPoint = 79;
							until ($breakPoint > length($geneSeq))
								{
								my $replaceString = "\n" . substr($geneSeq, $breakPoint, 1);
								substr($geneSeq, $breakPoint, 1) = $replaceString;
								$breakPoint += 80;
								}

						# Append sequence to file
							$geneStatus{$gene} = "stripped";
							$sequenceCount{$gene}++;
							$sequenceCount{$species{$accNum}}++;
							$sequenceCount{"all"}++;
							my $fastaFile = $gene."_gbs.fasta.txt";
							
							my $IOicon = ">";
								$IOicon .= ">" if ($sequenceCount{$gene} > 1);

							unless ($debug)


								{
								open (GENE, $IOicon, $fastaFile) or die "Cannot open file $fastaFile to write sequence to.";
									print GENE ">$species{$accNum} $accNum\n";
									print GENE "$geneSeq\n";
								close GENE;
								}

							# Update various statistical counters
								unless ($globalGenePresent{$gene})	# Gene counter
									{
									push @globalGeneList, $gene;
									$globalGenePresent{$gene} = 1;
									}

								unless ($speciesPresent{$species{$accNum}})	# Species counter
									{
									push @speciesList, $species{$accNum};
									$speciesPresent{$species{$accNum}} = 1;
									}

								$speciesGenePresent{$gene}{$species{$accNum}}++;	# Species-gene counter
									$speciesCount{$gene}++ if ($speciesGenePresent{$gene}{$species{$accNum}} == 1);
						
						# Print out summary stats for gene	
							if ($debug)
								{
								print "\t\t\tGene: $gene\n";
								printf "\t\t\t\tGene boundaries (%s bp)", length($geneSeq) - ($geneSeq =~ tr/\n//);
								print " (complemented)" if ($compStatus{$gene});
								print ":\n";
									for (my $i = 0; $i < scalar(@ { $geneStart{$gene} }); $i++)
										{
										print "\t\t\t\t\t$geneStart{$gene}[$i]\t$geneStop{$gene}[$i]\n";
										}
								print "\t\t\t\t\t$geneSeq\n";
								}
						}
					}
			}
close DATA;

# Print out summary stats
	print "\n\tSummary stats:\n";
	printf "\t\tTotal number of accessions processed: %s\n", scalar(@allAccNum);
	printf "\t\tTotal number of unique species: %s\n", scalar(@speciesList);
	printf "\t\tTotal number of unique sequences: %s\n", $sequenceCount{"all"};
	printf "\t\tTotal number of unique genes: %s\n", scalar(@globalGeneList);
	
	my $stripTime = time - $stripZero;
		print "\n\t\tTime taken to process file: $stripTime seconds\n";
		
	if (not @globalGeneList)
		{
		print  "\nNOTE: No genes stripped from file $gbFile; exiting\n";
		exit(0);
		}

# Reprocess results to 1) recheck whether stripped genes actually fall below user-set thresholds (geneCount procedure is a rough maximum), 2) pare down to required number of sequences per species, and 3) produce appropriate output files
	print "\nPostprocessing results ...\n";
		my $postZero = time;
		
	@globalGeneList = sort(@globalGeneList);
		$stripCount = scalar(@globalGeneList);
		
	foreach my $gene (@globalGeneList)
		{
		next if ($debug);
		print "\tProcessing gene $gene ($sequenceCount{$gene} sequences for $speciesCount{$gene} species)\n" if ($verbose);
		if ($sequenceCount{$gene} < $seqThreshold or $speciesCount{$gene} < $speciesThreshold)	# Gene does not meet threshold requirements; reject
			{
			$geneStatus{$gene} = "rejected";
				push @rejectList, $gene;
				$stripCount--;
			print "\t\tRejected; does not meet threshold requirements\n" if ($verbose);
			
			# Remove associated file
				my $rmString = $gene."_gbs.fasta.txt";
					$rmString =~ s/\s/\\ /g;
					$rmString =~ s/\(/\\(/g;
					$rmString =~ s/\)/\\)/g;
					$rmString =~ s/\'/\\'/g;
				system ("rm $rmString");
			}
		else
			{
			my (%infLength, %lengthList, %criticalLength);
			undef @speciesList;
			$geneStatus{$gene} = "stripped";
			if ($verbose)
				{
				print "\t\tStripped";
				print "; meets threshold requirements" if ($seqThreshold or $speciesThreshold);
				print "\n";
				}
			
			# Reload sequences from disk
				my $inputFile = $gene."_gbs.fasta.txt";
				setLineBreak($inputFile);
				open (FASTA, "<$inputFile") or die "Cannot open file $inputFile to read data for gene $gene\n";
					my (@accList, %speciesName, %fastaSeq, $species, $fastaAcc);

					undef %speciesPresent;
					while (<FASTA>)
						{
						chomp;
						if ($_ =~ "^>")
							{
 							($species, $fastaAcc) = split;
								$species =~ s/^>//;
								$fastaAcc =~ s/^\(//;
								$fastaAcc =~ s/\)$//;


							push @accList, $fastaAcc;
								$speciesName{$fastaAcc} = $species;
								$speciesPresent{$species}++;
								push @speciesList, $species if ($speciesPresent{$species} == 1);
							}
						else
							{
							$fastaSeq{$fastaAcc} .= $_;
							}
						}
				close FASTA;
				@accList = sort { $speciesName{$a} cmp $speciesName{$b} } keys %speciesName;
			
			# Pare down to desired number of sequences as needed
				if ($keepGene)
					{
					undef %lengthList;
					print "\t\tParing down to $keepGene best lengths for each of $speciesCount{$gene} species ...\n" if ($verbose);
					# Get informative length for each sequence
						foreach my $entry (@accList)
							{
							$infLength{$entry} = ($fastaSeq{$entry} =~ tr/ACGT//);
							push @ { $lengthList{$speciesName{$entry}} }, $infLength{$entry};
							}

					# Determine critical length and correct sequence count for number of deleted species
						foreach my $species (@speciesList)
							{
							print "\t\t\tSpecies: $species ($speciesPresent{$species} sequences)\n" if ($debug);
							next if ($speciesPresent{$species} <= $keepGene);	# Only process species for which gene has been sampled more than threshold
							@ { $lengthList{$species} } = sort {$b <=> $a} @ { $lengthList{$species} };	# Sort in descending order and get critical length
								$criticalLength{$species} = $lengthList{$species}[$keepGene-1];
	
							# Correct sequence count
								my $i = $keepGene;
									$i++ until ($i >= $speciesPresent{$species} or $lengthList{$species}[$i] < $criticalLength{$species});
									$sequenceCount{$gene} -= (scalar(@ { $lengthList{$species} }) - $i);
									
							if ($debug)
								{
								print "\t\t\t\tCritical length: $criticalLength{$species}\n";
								printf "\t\t\t\tNumber of sequences deleted: %s\n", scalar(@ { $lengthList{$species} }) - $i;
								print "\t\t\t\tTotal number of sequences for gene remaining: $sequenceCount{$gene}\n";
								}
							}
							
					print "\t\t\tFinal count: $sequenceCount{$gene} sequences\n" if ($verbose);
					
					# Recheck if paring has dropped sequence below sequence threshold
						if ($sequenceCount{$gene} < $seqThreshold)
							{
							print "\t\t\tNOTE: Gene ($sequenceCount{$gene} sequences) now falls below threshold of $seqThreshold sequences; no file created\n" if ($verbose);
							$geneStatus{$gene} = "rejected";
								push @rejectList, $gene;
								$stripCount--;

							# Remove associated file
								my $rmString = $gene."_gbs.fasta.txt";
									$rmString =~ s/\s/\\ /g;
									$rmString =~ s/\(/\\(/g;
									$rmString =~ s/\)/\\)/g;
									$rmString =~ s/\'/\\'/g;
								system ("rm $rmString");
							
							next;
							}
					}
			
			# Produce appropriate output files
				print "\t\tSaving output ...\n" if ($verbose);
				my %printCount;
				
				# Open files as needed
					my $sealFile = $gene."_gbs.seal";
						open (SEAL, ">$sealFile") or die "Cannot open Se-Al formatted file $sealFile for writing\n" if ($sealPrint);
					my $nexusFile = $gene."_gbs.nex";
						open (NEX, ">$nexusFile") or die "Cannot open nexus-formatted file $nexusFile for writing\n" if ($nexusPrint);
#					my $phytabFile = $gene."_gbs.phytab";
#						open (PHYTAB, ">$phytabFile") or die "Cannot open nexus-formatted file $phytabFile for writing\n" if ($phytabPrint);
#Currently assumes that only 1 file will be written for ease of Galaxy implementation
					my $phytabFile = "osiris_gbs.phytab";
						open (PHYTAB, ">$phytabFile") or die "Cannot open nexus-formatted file $phytabFile for writing\n" if ($phytabPrint);
					my $fastaFile = $gene."_gbs.fasta.txt";
						$fastaFile = $gene."_gbs.new.fasta.txt" if ($debug);
						open (FASTA, ">$fastaFile") or die "Cannot open fasta-formatted file $fastaFile for writing\n";
					
				# Print headers
					if ($sealPrint)
						{
						print SEAL "Database={\n";
						print SEAL "\tID='MLst';\n";
						print SEAL "\tOwner=null;\n";
						print SEAL "\tName=null;\n";
						print SEAL "\tDescription=null;\n";
						print SEAL "\tFlags=0;\n";
						print SEAL "\tCount=2;\n";
						print SEAL "\t{\n\t\t{\n";
						
						print SEAL "\t\t\tID='PAli';\n";
						print SEAL "\t\t\tOwner=1;\n";
						printf SEAL "\t\t\tName=\"%sgbs\";\n", $gene;
						print SEAL "\t\t\tDescription=null;\n";
						print SEAL "\t\t\tFlags=0;\n";
						print SEAL "\t\t\tNumSites=$maxLength{$gene};\n";
						print SEAL "\t\t\tType=\"Nucleotide\";\n";
						print SEAL "\t\t\tFeatures=null;\n";
						print SEAL "\t\t\tColourMode=1;\n";
						print SEAL "\t\t\tLabelMode=0;\n";
						print SEAL "\t\t\ttriplets=false;\n";
						print SEAL "\t\t\tinverse=true;\n";
						printf SEAL "\t\t\tCount=%s;\n", $sequenceCount{$gene};
						print SEAL "\t\t\t{\n";
						}

					if ($nexusPrint)
						{
						print NEX "#NEXUS\n\n";
						print NEX "Begin data;\n";
						printf NEX "\tDimensions ntax=%s nchar=%s;\n", $sequenceCount{$gene}, $maxLength{$gene};
						print NEX "\tFormat datatype=nucleotide gap=-;\n\n";
						print NEX "\tmatrix\n\n";
						}
						
				# Print sequence data
					my $i = 0;
					foreach my $entry (@accList)
						{
						next if (defined $criticalLength{$speciesName{$entry}} and $infLength{$entry} < $criticalLength{$speciesName{$entry}});
						
						$printCount{$speciesName{$entry}}++;
							$speciesName{$entry} .= "_$printCount{$speciesName{$entry}}" if ($printCount{$speciesName{$entry}} > 1);

						if ($sealPrint)
							{
							$i++;
							print SEAL "\t\t\t\t{\n";
							print SEAL "\t\t\t\t\tID='PSeq';\n";
							print SEAL "\t\t\t\t\tOwner=2;\n";
							print SEAL "\t\t\t\t\tName=\"$speciesName{$entry}\";\n";
							print SEAL "\t\t\t\t\tDescription=null;\n";
							print SEAL "\t\t\t\t\tFlags=0;\n";
							print SEAL "\t\t\t\t\tAccession=\"$entry\";\n";
							print SEAL "\t\t\t\t\tType=\"DNA\";\n";
							printf SEAL "\t\t\t\t\tLength=%s;\n", length($fastaSeq{$entry});
							print SEAL "\t\t\t\t\tSequence=\"$fastaSeq{$entry}\";\n";
							print SEAL "\t\t\t\t\tGeneticCode=1;\n";
							print SEAL "\t\t\t\t\tCodeTable=null;\n";
							print SEAL "\t\t\t\t\tFrame=1;\n";
							print SEAL "\t\t\t\t\tFeatures=null;\n";
							print SEAL "\t\t\t\t\tParent=null;\n";
							print SEAL "\t\t\t\t\tComplemented=false;\n";
							print SEAL "\t\t\t\t\tReversed=false;\n";
							print SEAL "\t\t\t\t}";
							print SEAL "," unless ($i == $sequenceCount{$gene});
							print SEAL "\n";
							}

						if ($nexusPrint)
							{
							my $gaps = "-" x ($maxLength{$gene} - length($fastaSeq{$entry}));
							print NEX "$speciesName{$entry}\t$fastaSeq{$entry}$gaps\n";
							}
						if($phytabPrint)
							{
							print PHYTAB "$speciesName{$entry} \t$gene \t $entry \t $fastaSeq{$entry} \n"
							}

						my $breakPoint = 79;
						until ($breakPoint > length($fastaSeq{$entry}))
							{
							my $replaceString = "\n" . substr($fastaSeq{$entry}, $breakPoint, 1);
							substr($fastaSeq{$entry}, $breakPoint, 1) = $replaceString;
							$breakPoint += 80;
							}
						if ($accessionFirst) #ADDED BY THO
							{
							print FASTA ">$entry", "\___", "$speciesName{$entry}\n";
							print FASTA "$fastaSeq{$entry}\n";
							}
						else
							{
							print FASTA ">$speciesName{$entry} ($entry)\n";
							print FASTA "$fastaSeq{$entry}\n";
							}
						}
				# Print footers
					if ($sealPrint)
						{
						print SEAL "\t\t\t};\n";
						print SEAL "\t\t},\n";
						print SEAL "\t\t{\n";
						print SEAL "\t\t\tID='MCoL';\n";
						print SEAL "\t\t\tOwner=1;\n";
						print SEAL "\t\t\tName=\"Genetic Codes\";\n";
						print SEAL "\t\t\tDescription=\"Custom Genetic Codes\";\n";
						print SEAL "\t\t\tFlags=0;\n";
						print SEAL "\t\t\tCount=0;\n";
						print SEAL "\t\t}\n";
						print SEAL "\t};\n";
						print SEAL "};\n";
						}
					if ($nexusPrint)
						{
						print NEX "\t;\nend;\n";
						}
				
				# Close files as appropriate
					close SEAL if ($sealPrint);
					close NEX if ($nexusPrint);
					close FASTA if ($keepGene);
			}
		}

	my $postTime = time - $postZero;
	print "\n" if ($verbose);
	print "\tTime taken: $postTime seconds\n";

# Print out final summary stats / files
	print "\nFinal summary statistics\n";
	
	# Print file of stripped genes
		open (STRIP, ">$stripFile") or die "Cannot write summary of stripped genes to $stripFile.\n";
			print STRIP "The following $stripCount genes were stripped successfully from file $gbFile\n";
			print STRIP "\nGene\tNo. of sequences\tNo. of species\n\n";
	
			foreach my $gene (@globalGeneList)
				{
				print STRIP "$gene\t$sequenceCount{$gene}\t$speciesCount{$gene}\n" unless ($geneStatus{$gene} eq "rejected");
				}
		close STRIP;
		
		print "\tSummaries of $stripCount genes that were stripped from $gbFile have been written to $stripFile\n";

	# Print file of rejected genes
		if (@rejectList)
			{
			@rejectList = sort @rejectList;
			open (REJECT, ">$rejectFile") or die "Cannot write summary of rejected genes to $rejectFile.\n";
				printf REJECT "The following %s genes do not meet user-defined threshold(s) of", scalar(@rejectList);
					print REJECT " $seqThreshold sequences" if ($seqThreshold);
					print REJECT " or" if ($seqThreshold and $speciesThreshold);
					print REJECT " $speciesThreshold species" if ($speciesThreshold);
					print REJECT "; an asterisk indicates that values given are rough maximums\n";
				print REJECT "\nGene\tNo. of sequences\tNo. of species\n\n";

				foreach my $gene (@rejectList)
					{
					print REJECT "$gene\t";
					if (defined $sequenceCount{$gene})	# Need to do alternative versions depending on whether gene was properly counted or not
						{
						print REJECT "$sequenceCount{$gene}\t";
						}
					else
						{
						print REJECT "$quickSequenceCount{$gene} \*\t";
						}
					if (defined $speciesCount{$gene})
						{
						print REJECT "$speciesCount{$gene}\n";
						}
					else
						{
						print REJECT "$quickSpeciesCount{$gene} \*\n";
						}
					}
			close REJECT;

			printf "\tSummaries of %s genes that did not meet threshold(s) have been written to $rejectFile\n", scalar(@rejectList);
			}

	if ($debug)
		{
		print "\n\tSummary stats for each processed gene:\n";
		foreach my $gene (@globalGeneList)
			{
			print "\t\tGene: $gene ($geneStatus{$gene})\n";
			print "\t\t\tTotal number of sequences: $sequenceCount{$gene}\n";
			print "\t\t\tTotal number of unique species: $speciesCount{$gene}\n";
			}
		}

exit(0);

# Subroutines used in the script

sub setLineBreak	# Check line breaks of input files and set input record separator accordingly
	{
	my $gbFile = shift;
	$/ ="\n";
	open (IN, "<$gbFile") or die "Cannot open $gbFile to check form of line breaks.\n";
		while (<IN>)
			{
			if ($_ =~ /\r\n/)
				{
				print "\tDOS line breaks detected ...\n" if ($debug);
				$/ ="\r\n";
				last;
				}
			elsif ($_ =~ /\r/)
				{
				print "\tMac line breaks detected ...\n" if ($debug);
				$/ ="\r";
				last;
				}
			else
				{
				print "\tUnix line breaks detected ...\n" if ($debug);
				$/ ="\n";
				last;
				}
			}
	close IN;
	}
	
sub complement	# Outputs complementary sequence to one provided
	{
	my $tempSeq = shift;

	my $compSeq;
		for (my $nt = 0; $nt < length($tempSeq); $nt++)
			{
			if (not defined $complement{substr($tempSeq, $nt, 1)})
				{
				$compSeq .= "?";
				}
			else
				{
				$compSeq .= $complement{substr($tempSeq, $nt, 1)};
				}
			}
        my $revCompSeq = reverse($compSeq);
	return $revCompSeq;
	}

sub geneSynonyms	# Define gene synonyms; originally compiled by Robin Beck
	{
	#Opsin gene added by THO
	$synonym{$_} = "opsin" foreach ("opsin", "anceropsin", "blop", "blue-sensitive_opsin", "blue-sensitive_opsin_precursor", "blue-sensitive_rhodopsin", "blue-sensitive_visual_pigment", "blue_opsin", "bluerh", "boceropsin", "buvops", "compound_eye_opsin_bcrh1", "compound_eye_opsin_bcrh2", "lateral_eye_opsin", "locus_opsin_1", "locust_opsin_2", "long-wavelenght_opsin", "long-wavelength_like_opsin", "long-wavelength_opsin", "long-wavelength_rhodopsin", "long-wavelength_sensitive_opsin_1", "long-wavelength_sensitive_opsin_2", "long_wave_opsin", "long_wavelength-sensitive_opsin", "long_wavelength-sensitive_rhodopsin", "long_wavelength-sensitive_visual_pigment", "long_wavelength_opsin", "long_wavelength_sensitive_opsin_1", "long_wavelength_sensitive_opsin_2", "lop2", "lw_opsin", 
"ocellar_opsin", "opsin", "opsin_2", "opsin_bcrh1", "opsin_bcrh2", 
"opsin_rh1", "opsin_rh3", "opsin_rh4", "piceropsin", "pteropsin", 
"rh1_opsin", "rh2_opsin", "rh3_opsin", "rh4_opsin", "rh6_rhodopsin", 
"rhodopsin", "rhodopsin_1", "rhodopsin_2_cg16740-pa", "rhodopsin_3", 
"rhodopsin_3_cg10888-pa", "rhodopsin_4", "rhodopsin_4_cg9668-pa", 
"rhodopsin_5", "rhodopsin_5_cg5279-pa", "rhodopsin_6", 
"rhodopsin_6_cg5192-pb", "rhodopsin_7_cg5638-pa", 
"rhodopsin_long-wavelength", "short_wavelength-sensitive_opsin", 
"ultraviolet-sensitive_opsin", "ultraviolet-sensitive_rhodopsin", 
"ultraviolet-sensitive_visual_pigment", "ultraviolet_sensitive_opsin", 
"uv-sensitive_opsin", "uv-wavelength_like_opsin", "uv_opsin", "uvop", 
"uvrh", "uvrh1", "amblop", "amuvop", 
"lwrh", "lwrh1", "lwrh2", "lw", "lw-rh", "lw_rh", "ninae", "ninae-pa", 
"op", "ops", "ops1", "opsin_1", "opsin_3", "rh", "rh1", "rh2", "rh3", 
"rh2-pa", "rh3-pa", "rh4", "rh4-pa", "rh5", "rh6", "rh7", "rho", 
"visual_pigment", "long-wavelength_rodopsin", 
"long_wavelength_rhodopsin", "short-wavelength_rhodopsin");

	#Chloroplast genes added by THO***
		$synonym{$_} = "rbcl" foreach ("rbcl", "large_subunit_of_riblose-15-bisphosphate_carboxylase-oxygenase","larger_subunit_of_rubisco", 
						"rubisco_large_subunit", "ribulose_15-bisphosphate_carboxylase_large_subunit", "ribulosebiphosphate_carboxylase_large_subunit"); 
		$synonym{$_} = "matk" foreach ("matk", "maturase_k", "maturase");

	# Mitochondrial genes
		$synonym{$_} = "mtatp6" foreach ("mtatp6", "atp synthase 6", "atp synthase a chain", "atp6", "mt-atp6",
										"atpase subunit 6", "atpase6", "atpase 6", "atpase6 protein", "atp6",
										"atp synthase f0 subunit 6", "atp synthase a chain protein 6",
										"atp synthase subunit 6", "atpase6", "atp6", "atp sythase subunit 6",
										"atpase subunit 6, atpase6", "f0-atpase subunit 6", "f1atpase subunit 6",
										"f1-atpase subunit 6", "atpase 6", "atpase 6 gene", "atpase subunit 6 (atp6)",
										"atpase subunit 6 (atpase6)");
		$synonym{$_} = "mtatp8" foreach ("mtatp8", "atp synthase 8", "atp synthase protein 8", "atpase subunit 8",
										"a6l", "atp8", "mt-atp8", "atpase subunit 8", "atpase8", "atpase 8",
										"atpase8 protein", "atp8", "atp synthase f0 subunit 8", "atp synthase protein 8",
										"atpase-8", "atp synthase subunit 8", "atpase8", "atp8", "atp sythase subunit 8",
										"atpase subunit8", "f0-atpase subunit 8", "f1 atpase subunit 8", "f1-atpase subunit 8",
										"protein a6l", "atpase 8", "atpase subunit 8 (atp8)", "atpase subunit 8 (atpase8)",
										"atpase 8 gene", "atpase subunit 8 (atp8)", "atpase subunit 8 (atpase8)");
		$synonym{$_} = "mtco1" foreach ("mtco1", "cytochrome c oxidase i", "coi", "mt-co1", "co i", "ccoi", "cox 1", "cox i",
										"coi", "cytochrome c oxidase subunit i", "cytochrome oxidase subunit i", "cox1",
										"cytochrome c oxidase subunit 1", "cytochrome oxidase subunit 1", "co1",
										"cytochrome c oxidase polypeptide i", "cytochrome oxidase i", "cox-1",
										"cytochrome oxidase c subunit 1", "cox1", "co i", "co1", "cox 1", "coxi",
										"cytochrome c oxidase polypeptide i", "cytochrome c-oxidase subunit 1",
										"cytochrome oxidase c subunit i", "cytochrome-c oxidase i", "cytochrome c1 mrna",
										"cytochrome c1", "cytochrome c oxidase subunit 1 (coxi)", "cytochrome c oxidase chain i",
										"cytochrome c1 (aa 1-241)", "cytochrome c oxidase polypeptide 1",
										"cytochrome c oxidase subunit 1 (coi)", "cytochrome c-1");
		$synonym{$_} = "mtco2" foreach ("mtco2", "cytochrome c oxidase ii", "cytochrome c oxidase polypeptide ii", "coii",
										"mt-co2", "coii", "cytochrome c oxidase subunit ii", "cytochrome oxidase subunit ii",
										"cytochrome oxidase subunit 2", "cytochrome c oxidase subunit 2",
										"cytochrome c oxidase ii", "cox2", "co2", "cytochrome c oxidase polypeptide ii",
										"cytochrome oxidase ii", "cox2", "cytochrome oxidase c subunit 2",
										"cytochrome-c oxidase", "co ii", "co2", "cox 2", "cytochrome c oxidase polypeptide ii",
										"cytochrome c-oxidase subunit 2", "cytochrome oxidase c subunit ii",
										"cytochrome oxidase subunit2", "cytochrome-c oxidase ii", "cytochrome c oxidase subunit 2 (coxii)",
										"cytochrome c oxidase subunit 2 (coii)", "cytochrome c oxidase chain ii",
										"cytochrome c oxidase polypeptide 2", "cytochrome c oxidase ii subunit",
										"cytochrome c oxidase ii mrna");
		$synonym{$_} = "mtco3" foreach ("mtco3", "cytochrome c oxidase iii", "coiii", "mt-co3", "co iii", "ccoiii", "cox3",
										"cox iii", "cytochrome oxidase subunit iii", "coiii", "cox iii",
										"cytochrome c oxidase subunit iii", "cytochrome oxidase subunit 3", "cox3",
										"cytochrome c oxidase subunit 3", "co3", "cytochrome c oxidase polypeptide iii",
										"cox3", "cytochrome oxidase c subunit 3", "cytochrome oxidase iii",
										"cytochrome c oxidase polypeptide iii", "co iii", "coiii", "coiii protein", "coxiii",
										"cytochrome c oxidase subunit iii, coiii", "cytochrome c-oxidase subunit 3",
										"cytochrome c-oxidase subunit three", "cytochrome oxidase c subunit iii",
										"cytochrome-c oxidase iii", "cytochrme c oxidase iii", "cytochrome c oxidase subunit 3 (coiii)",
										"cytochrome oxidase subunit iii type 2");
		$synonym{$_} = "mtnd1" foreach ("mtnd1", "nd1", "nadh dehydrogenase 1", "nadh dehydrogenase, subunit 1 (complex i)",
										"mt-nd1", "nadh-ubiquinone oxidoreductase chain 1", "nd1", "nadh1", "nd1",
										"nadh-ubiquinone oxidoreductase chain 1", "nadh-ubiquinone oxidoreductase subunit 1",
										"nadh dehydrogenase i", "nadh dehydrogenase 1", "nadh subunit 1",
										"nadh dehydrogenase subnuit 1", "nadh1", "nadh1 protein",
										"nadh-ubiquinone oxidoreductase subunit i", "nadh dehydrogenase subunit 1",
										"nadh dehydrogenase subunit 1 (nd1)", "nadh dehydrogenase subunit 1 type 2",
										"nadh dehydrogenase subunit 1 type 1", "nadh dehydrogenase subunit i", "nadh dehydrogenase i",
										"nadh dehydrogenase subnit 1", "nadh dehydrogenase ubiquinone 1 alpha",
										"nadh dehydrogenase subunit i");
		$synonym{$_} = "mtnd2" foreach ("ndh-u1", "mtnd2", "nadh dehydrogenase 2", "nadh dehydrogenase, subunit 2 (complex i)",
										"nadh-ubiquinone oxidoreductase chain 2", "nd2", "mt-nd2", "urf2",
										"nadh dehydrogenase subunit 2", "nd2", "nadh2",
										"nadh-ubiquinone oxidoreductase chain 2", "nadh dehydrogenase 2", "nadh subunit 2",
										"nadh-ubiquinone oxidoreductase subunit 2", "nadh dehydrogenase subnuit 2",
										"nadh deydrogenase subunit 2", "nadh2", "nadh-ubiquinone oxidoreductase subunit ii",
										"nd2", "nadh2 protein", "nadh dehydrogenase subunit 2 (nd2)", "nadh dehydrogenase subunit ii",
										"nadh dehydrogenase ii", "nadh dehydrogense 2", "nadh dehydrogenase subunit ii");
		$synonym{$_} = "mtnd3" foreach ("mtnd3", "ndh-u3", "nadh dehydrogenase 3", "nd3", "nadh3",
										"nadh-ubiquinone oxidoreductase chain 3",
										"nadh-ubiquinone/plastoquinone oxidoreductase chain 3", "mt-md3", "urf3",
										"nadh dehydrogenase subunit 3", "nd3", "nadh3",
										"nadh-ubiquinone oxidoreductase chain 3", "nadh subunit 3", "nadh3", "nadh3 protein",
										"nadh-ubiquinone oxidoreductase subunit 3", "nadh dehydrogenase subnuit 3",
										"nadh-ubiquinone oxidoreductase subunit iii", "nadh3 protein",
										"nadh dehydrogenase subunit 3 (nd3)", "nadh dehydrogenase subunit iii", "nadh dehydrogenase iii",
										"nadh dehydrogenase subunit iii");
		$synonym{$_} = "mtnd4" foreach ("mtnd4", "ndh-u4", "nadh dehydrogenase 4", "nadh-ubiquinone oxidoreductase chain 4", "mt-nd4",
										"nd4", "urf4", "nadh:ubiquinone oxidoreductase subunit 4 (chain m)",
										"nadh dehydrogenase subunit 4", "nd4", "nadh4",
										"nadh-ubiquinone oxidoreductase chain 4", "nadh subunit 4", "nadh4",
										"nadh-ubiquinone oxidoreductase subunit 4", "nadh dehydrogenase subunit4",
										"nadh-ubiquinone oxidoreductase subunit iv", "nadh4 protein",
										"nadh dehydrogenase subunit 4 (nd4)", "nadh dehydrogenase subunit iv", "nadh dehydrogenase iv",
										"nadh dehydrogenase subunit iv");
		$synonym{$_} = "mtnd4l" foreach ("mtnd4l", "ndh-u4l", "nadh dehydrogenase 4l", "nadh dehydrogenase, subunit 4l (complex i)",
										"nadh-ubiquinone oxidoreductase chain 4l", "nd4l", "mt-md4l", "urf4l",
										"nadh dehydrogenase subunit 4l", "nd4l", "nadh4l", "nadh-ubiquinone oxidoreductase chain 4l",
										"nadh subunit 4l", "nadh4l", "nadh-ubiquinone oxidoreductase subunit 4l",
										"nadh dehydrogenase subunit 4 l", "nadh4l protein", "nadh dehydrogenase subunit 4l (nd4l)",
										"nadh dehydrogenase subunit ivl", "nadh dehydrogenase ivl", "nadh dehydrogenase subunit ivl");
		$synonym{$_} = "mtnd5" foreach ("mtnd5", "ndh-u5", "nadh dehydrogenase 5", "nd5", "nadh-ubiquinone oxidoreductase chain 5",
										"mt-nd5", "urf5", "nadh dehydrogenase subunit 5", "nadh5", "nd5",
										"nadh dehydrogenase-5", "nadh-ubiquinone oxidoreductase chain 5",
										"nadh dehydrogenase 5", "nadh 5", "nadh subunit 5",
										"nadh-ubiquinone oxidoreductase subunit 5", "nadh5", "nadh-dehydrogenase subunit 5",
										"nadh-dehydrogenase subunit 5, nd5", "nadh-ubiquinone oxidoreductase subunit v",
										"nadh5 protein", "nadh dehydrogenase subunit 5 (nd5)", "nadh dehydrogenase subunit v",
										"nadh dehydrogenase v", "nadh dehydrogenase subunit v");
		$synonym{$_} = "mtnd6" foreach ("mtnd6", "ndh-u6", "nadh dehydrogenase 6", "nd6", "nadh6",
										"nadh-ubiquinone oxidoreductase chain 6",
										"nadh-ubiquinone/plastoquinone oxidoreductase chain 6", "mt-md6", "urf6",
										"nadh dehydrogenase subunit 6", "nd6", "nadh6",
										"nadh-ubiquinone oxidoreductase chain 6", "nadh subunit 6",
										"nadh-ubiquinone oxidoreductase subunit 6", "nadh dehydrogenase 6", "nadh6",
										"nadh-ubiquinone oxidoreductase subunit vi", "nadh6 protein",
										"nadh dehydrogenase subunit 6 (nd6)", "nadh dehydrogenase subunit vi", "nadh dehydrogenase vi",
										"nadh dehydrogenase subunit vi");
		$synonym{$_} = "mtrnr1" foreach ("mtrnr1", "mt-rnr1", "12s rna", "12s rrna", "12s ribosomal rna", "12s rrna",
										"small subunit ribosomal rna", "12 rrna", "12 s ribosomal rna", "12s rrna, mtssu rrna",
										
"12s rrna gene", "12s small subunit ribosomal rna", 
"mitochondrial ssu ribosomal rna", "rrns",
										"ssu ribosomal rna", "12s","12s_ribosmal_rna", "s-rrna", "srrna", "ssu");
		$synonym{$_} = "mtrnr2" foreach ("16s", "16s_ribosomal_rna","l-rrna","lrrna","lsu", "mtrnr2", "mt-rnr2", "16s rrna", "16s rna", "16srrna", "16s ribosomal rna", "16s rrna", "rrna large subunit",
										"large subunit ribosomal rna", "16 s ribosomal rna", "16s mitochondrial ribosomal rna",
										
"mitochondrial 16s ribosomal rna", "16s large subunit ribosomal 
rna", "16s_ribosmal_rna","rrnl");

	# Nuclear genes
		$synonym{$_} = "adora3" foreach ("adora3", "adenosine a3 receptor", "a3ar", "ara3", "gpcr2", "adora3", "adenosine a3 receptor",
										"a3 adenosine receptor", "adenosine-3 receptor");
		$synonym{$_} = "adra2b" foreach ("adra2b", "adrenergic, alpha-2b-, receptor", "adra2l1", "adrarl1", "adra2rl1",
										"alpha-2b adrenergic receptor", "alpha-2b adrenoceptor", "subtype c2", "[a]2b", "adra-2b",
										"alpha2-c2", "alpha2b", "subtype alpha2-c2", "alpha-2-adrenergic receptor-like 1",
										"alpha adrenergic receptor 2b", "adra2b", "alpha 2b adrenergic receptor", "adra2b",
										"alpha adrenergic receptor subtype 2b", "alpha-2b-adrenergic receptor",
										"alpha adrenergic receptor, subtype 2b", "alpha-2b adrenergic receptor", "alpha-2b adrenoceptor");
		$synonym{$_} = "adrb2" foreach ("adrb2", "adrenergic, beta-2-, receptor, surface", "adrb2r", "adrbr", "beta-2 adrenergic receptor",
										"b2ar", "bar", "beta-2 adrenoceptor", "catecholamine receptor", "adrb-2", "badm", "beta 2-ar",
										"gpcr7", "adrb2", "beta-2 adrenergic receptor", "beta-2-adrenergic receptor",
										"adrenergic receptor beta 2", "beta 2 adrenergic receptor", "beta2-adrenergic receptor",
										"beta2 adrenergic receptor");
		$synonym{$_} = "apob" foreach ("apob", "apolipoprotein b", "ag(x) antigen", "apolipoprotein b-100 precursor", "apo b-100",
										"apolipoproteinb-48", "apo b-48", "fldb", "apob-100", "apolipoprotein b", "apolipoprotein b 100",
										"apob", "apolipoprotein b-100", "apob");
		$synonym{$_} = "app" foreach ("app", "amyloid beta (a4) precursor protein", "protease nexin-ii, alzheimer disease", "ad1",
										"human mrna for amyloid a4 precursor of alzheimer's disease", "appi", "beta-amyloid protein",
										"beta-app", "a-beta", "a4", "cvap", "aaa", "abeta", "amyloid beta-peptide",
										"amyloid beta (a4) precursor protein", "adap", "appican", "betaapp", "app",
										"amyloid beta precursor protein", "amyloid precursor protein", "amyloid beta precursor b1",
										"beta amyloid protein precursor", "beta-a4");
		$synonym{$_} = "atp7a" foreach ("atp7a", "atpase, cu++ transporting, alpha polypeptide (menkes syndrome)", "mnk",
										"copper-transporting atpase 1", "copper pump 1", "menkes disease-associated protein", "mc1",
										"copper binding p-type atpase 7a", "blo", "blotchy", "br", "brindled", "menkes protein", "mo",
										"mottled", "mk", "ohs", "atp7a", "copper transporting atpase",
										"atpase, cu++ transporting, alpha polypeptide", "menkes syndrome protein");
		$synonym{$_} = "bdnf" foreach ("bdnf", "brain-derived neurotrophic factor", "brain-derived neurotrophic factor precursor", "bdnf",
										"brain-derived neurotrophic factor", "brain derived neurotrophic factor",
										"brain-derived neurotrophic factor mature peptide", "bdnf");
		$synonym{$_} = "bmi1" foreach ("bmi1", "b lymphoma mo-mlv insertion region (mouse)",
										"murine leukemia viral (bmi-1) oncogene homolog", "polycomb complex protein bmi-1", "rnf51",
										"mgc12685", "oncogene bmi-1", "bmi1", "bmi-1", "bmi-1 protein", "oncoprotein bmi-1");
		$synonym{$_} = "brca1" foreach ("brca1", "breast cancer 1, early onset", "pscp", "papillary serous carcinoma of the peritoneum",
										"breast cancer type 1 susceptibility protein", "breast and ovarian cancer susceptibility gene",
										"rnf53", "breast-ovarian cancer, included", "brca1", "brca1",
										"breast and ovarian cancer susceptibility protein", "brca1 protein",
										"breast and ovarian cancer susceptibility");
		$synonym{$_} = "chrna1" foreach ("chrna1", "cholinergic receptor, nicotinic, alpha polypeptide 1 (muscle)", "chrna",
										"acetylcholine receptor protein, alpha chain precursor", "achra", "achr-1", "acra", "chrna1",
										"nicotinic cholinergic receptor alpha polypeptide", "chrna1", "achr", "achr prepeptide", "chrna",
										"neuronal nicotinic acetylcholine receptor alpha", "nicotinic acetylcholine recepter alpha-subunit");
		$synonym{$_} = "cftr" foreach ("cftr", "cystic fibrosis transmembrane conductance regulator, atp-binding cassette (sub-family c,
										member 7)", "mrp7", "cf", "abcc7", "abc35", "cftr", "cystic fibrosis transmembrane conductance",
										"cftr chloride channel");
		$synonym{$_} = "cnr1" foreach ("cnr1", "cannabinoid receptor 1 (brain)", "cnr", "cb1", "cb-r", "cann6", "cb>1<", "cb1k5", "cb1a",
										"central cannabinoid receptor", "cnr1", "cannabinoid receptor 1", "cb1", "cb1 cannabinoid receptor",
										"cb1 cannabinoid receptor", "cbr");
		$synonym{$_} = "crem" foreach ("crem", "camp responsive modulator", "crea", "camp-responsive element modulator, alpha isoform",
										"crem", "camp responsive element moderator", "camp responsive element modulator",
										"camp-responsive element moderator");
		$synonym{$_} = "edg1" foreach ("edg1", "endothelial differentiation, sphingolipid g-protein-coupled receptor, 1", "d1s3362", "edg-1",
										"1pb1", "s1p1", "ecgf1", "chedg1", "sphingosine 1-phosphate receptor edg1",
										"g protein-coupled sphingolipid receptor", "edg1");
		$synonym{$_} = "ghr" foreach ("ghr", "growth hormone receptor", "growth hormone receptor precursor", "gh receptor",
										"serum binding protein", "growth hormone receptor", "ghr", "growth hormone receptor precursor",
										"ghr", "bovine growth hormone receptor", "mature growth hormone receptor");
		$synonym{$_} = "pgk1" foreach ("pgk1", "phosphoglycerate kinase 1", "pgk-1", "primer recognition protein 2", "prp 2", "pgka",
										"pgk-1", "pgk-1", "phosphoglycerate kinase 1");
		$synonym{$_} = "plcb4" foreach ("plcb4", "phospholipase c, beta 4",
										"1-phosphatidylinositol-4,5-bisphosphate phosphodiesterase beta 4", "plc-beta-4",
										"phospholipase c-beta-4", "plcb4", "phospholipase c beta 4");
		$synonym{$_} = "pnoc" foreach ("pnoc", "prepronociceptin", "propronociceptin", "nociceptin precursor", "orphanin fq", "ppnoc", "ofq",
										"n/ofq", "n23k", "npnc1", "ofq/n", "proorphanin", "pnoc", "prepronociceptin",
										"nociceptin/orphanin fq precursor");
		$synonym{$_} = "prkci" foreach ("prkci", "protein kinase c, iota", "pkci", "dxs1179e", "npkc-iota", "protein kinase c iota",
										"prkci");
		$synonym{$_} = "prnp" foreach ("prnp", "(prion protein (p27-30) (creutzfeld-jakob disease, gerstmann-strausler-scheinker syndrome, fatal familial insomnia)",
										"cjd", "major prion protein precursor", "prp", "prp27-30", "prp33-35c", "ascr", "prip", "gss",
										"prn-i", "prn-p", "prpc", "prpsc", "sinc", "mgc26679", "cd230 antigen", "prion-related protein",
										"prion protein", "prp", "prnp", "prion protein precursor", "prion protein", "prnp", "prp",
										"prion protein prp", "prion protein precursor prp", "prp", "prp", "prp gene", "greater kudu prp",
										"major prion protein", "prion protein variant 110p", "prion protein variant 143r",
										"prion protein variant 240s", "prion protein variant 37v", "prion protein variant 37v/240s",
										"prion protein, prp", "prnp", "prmp");
		$synonym{$_} = "rag1" foreach ("rag1", "recombination activating gene 1", "v(d)j recombination activating protein 1", "rag-1",
										"recombination activating protein 1", "rag1", "rag-1", "recombination activating gene-1",
										"recombination activating gene 1", "recombination-activating gene 1", "rag1", "rag 1", "rag1");
		$synonym{$_} = "rag2" foreach ("rag2", "recombination activator protein 2", "recombination activating protein 2", "rag-2", "rag 2",
										"recombination activating protein", "recombination activating gene-2", "rag2",
										"recombination activating protein 2", "rag2", "rag2", "rag2", "rag-2", "rag-2 protein",
										"recombinase activating gene 2", "recombination activating protein 2");
		$synonym{$_} = "rbp3" foreach ("rbp3", "retinol-binding protein 3, interstitial", "interphotoreceptor retinoid-binding protein precursor",
										"irbp", "interstitial retinol-binding protein", "rbp-3", "interphotoreceptor retinoid binding protein",
										"irbp", "irbp", "interphotoreceptor retinoid-binding protein", "interphotorecepter retinoid binding protein",
										"interphotoreceptor retinoid binding protein (irbp)", "irbp mrna");
		$synonym{$_} = "tnf" foreach ("tnf", "tumor necrosis factor (tnf superfamily, member 2)", "tnfa", "dif", "tnfsf2",
										"tumor necrosis factor (cachectin)", "tumor necrosis factor, alpha (cachectin)",
										"tumor necrosis factor", "tnf-alpha", "cachectin", "tnfsf1a", "apc1 protein",
										"tnf, monocyte-derived", "tnf, macrophage-derived", "tnf superfamily, member 2",
										"tumor necrosis factor-alpha", "tumor necrosis factor alpha", "tnfa", "tnf-alpha",
										"tumor necrosis factor alpha precursor", "tnfa", "tnfa", "tnfalpha", "tumour necrosis factor alpha",
										"tnf alpha", "tnf-a", "tumor necrosis factor alpha, tnf-alpha", "bovine tumor necrosis factor alpha",
										"tumor necrosis factor alpha (cachetin)", "tumor necrosis factor-alpha precursor");
		$synonym{$_} = "tp53" foreach ("tp53", "tumor protein p53 (li-fraumeni syndrome)", "p53", "cellular tumor antigen p53",
										"phosphoprotein p53", "trp53", "transformation related protein 53", "p53", "p53 protein", "p53",
										"tp53", "tumor suppressor p53", "53 kda phosphoprotein", "insulin recptor substrate p53 short form",
										"p53 gene product", "p53 tumor suppressor gene", "p53 tumor suppressor protein",
										"tumor suppressor p53 phosphoprotein");
		$synonym{$_} = "ttr" foreach ("ttr", "transthyretin (prealbumin, amyloidosis type i)", "palb", "tthy", "tbpa", "attr",
										"transthyretin", "transthyretin precursor", "transthyretin subunit", "ttr", "ttr");
		$synonym{$_} = "tyr" foreach ("tyr", "tyrosinase (oculocutaneous albinism ia)", "ocaia", "monophenol monooxygenase",
										"tumor rejection antigen ab", "sk29-ab", "lb24-ab", "albino", "c", "skc35", "oca1a", "tyrosinase",
										"tyr", "tyrosinase precursor", "truncated tyrosinase", "tyr", "tyr");
		$synonym{$_} = "vwf" foreach ("vwf", "von willebrand factor", "f8vwf", "coagulation factor viii", "von willebrand factor", "vwf",
										"von willebrand factor", "vwf", "vwf", "vwf", "von willebrand factor, vwf", "von willebrand factor precursor",
										"von willebrand factor precursor", "wf");
		$synonym{$_} = "zfx" foreach ("zfx", "zinc finger protein, x-linked", "zinc finger x-chromosomal protein", "zfx", "zinc finger protein zfx",
										"zfx", "zinc finger protein zfx", "x-linked zinc finger protein zfx", "zfx", "x-linked zinc finger protein",
										"zinc finger x-linked protein 1", "zinc finger x-linked protein 2", "zinc finger protein x linked",
										"zfx1", "zfx-1", "zfx2", "zfx-2","zfx protein mrna", "zinc finger protein zfx isoform 4",
										"zfx product, isoform 1", "zfx product, isoform 2", "zfx product, isoform 3",
										"zfx product, isoform 4", "zfx protein");
		$synonym{$_} = "zfy" foreach ("zfy", "zinc finger protein, y-linked", "zinc finger y-chromosomal protein", "zfy", "zfy",
										"zinc finger protein zfy", "zinc finger protein zfy", "y-linked zinc finger protein",
										"y-linked zinc finger protein zfy", "zinc finger protein y linked",
										"y-linked zinc finger protein 2", "y-linked zinc finger protein 1", "zfy1", "zfy-1", "zfy2",
										"zfy-2");

	# Additional genes
		$synonym{$_} = "18s_rrna" foreach ("18s rrna", "18s_rrna", "18s ribosomal rna", 
										"18 rrna", "18 s ribosomal rna", "18 s ribosomal rna", "18s small subunit ribosomal rna", 
										"18s_ribsomal_rna", "18s_ribosomal_rna_a-type", "18s_rna_gene", "18s_rrna", "nuclear_18s_ribosomal_rna");

		$synonym{$_} = "28s_rrna" foreach ("28s_large_subunit_ribosomal_rna", "28s rrna", "28s_rrna", "28s ribosomal rna", "28 rrna", "put. 28S ribosomal RNA", "28 s ribosomal rna", "28 s ribosomal rna",
										"rrna-28s_rrna-related", "28s ribosomal rna v region");
		$synonym{$_} = "5.8s" foreach ("5.8s_ribosomal_rna", "5.8s_rrna", "5s_ribosomal_rna","5s_rrna","5.8s_rna_gene");    
		$synonym{$_} = "mtcyb" foreach ("cob", "mtcyb", "cyt b", "cytb", "cytochrome b", "cytochrome b", "cytb", "cyt b", "cytb", "cytb",
										"'cytochrome b'", "cytochrome-b", "skeletal muscle cytochrome b", "cyt-b", "cyt.b", "cyto b",
										"cytob", "cytb1", "cytb2", "cytb3", "cytb4", "cytochrome b light strand", "cyt-b", "cytob",
										"cyt.b", "cytb", "cytochrome b", "mitochondrial cytochrome b", "cytochrome b protein",
										"duodenal cytochrome b");
		$synonym{$_} = "alpha_lactalbumin" foreach ("alpha lactalbumin", "alpha_lactalbumin", "alpha-lactalbumin", "alpla lactalbumin",
										"alpah lactalbumin", "a-lacta", "a-lacta", "lactalbumin, alpha-", "lactalbumin, alpha",
										"mature alpha-lactalbumin");
		$synonym{$_} = "alpha-s1-casein" foreach ("alpha-s1-casein", "as1-casein", "alpha s1 casein", "alpha s1-casein b",
										"alpha s1 casein", "alpha s1-casein", "alpha(s1) casein", "alpha-s1 casein", "alpha-s1 casein",
										"alpha-s1-casein", "alpha-s1-casein mrna", "alphas1-casein", "as1-casein", "alfas1-casein",
										"alpha-sl casein", "casein alpha s1");
		$synonym{$_} = "alpha-s2-casein" foreach ("alpha-s2-casein", "as2-casein", "alpha s2 casein", "alpha s2-casein b",
										"alpha s2 casein", "alpha s2-casein", "alpha(s2) casein", "alpha-s2 casein", "alpha-s2 casein",
										"alpha-s2-casein", "alpha-s2-casein mrna", "alphas2-casein", "as2-casein",
										"alpha s2a casein (aa 1 to 165)", "alpha s2a casein (aa 1 to 167)", "alph as2-casein",
										"alpha(s2)-casein");
		$synonym{$_} = "b-casein" foreach ("b-casein", "beta casein", "beta-casein", "beta-casein (aa 1 - 213)", "beta-casein a3",
										"beta-casein precursor", "beta casein precursor", "beta-casein variant a2",
										"beta-casein variant i", "mat. beta-casein", "casein beta", "casein b", "casein b (aa 1-183) pgpk48",
										"mature beta-casein (aa 1 to 226)");
		$synonym{$_} = "k-casein" foreach ("k-casein", "kappa casein", "kappa-casein", "kappa-casein precursor", "kappa casein precursor",
										"casein kappa", "kappa-cas", "kappa-casein long form", "kappa-casein mature peptide",
										"kappa-casein short form");
		$synonym{$_} = "sry" foreach ("sry", "sex determining factor sry", "sex determining region y protein", "sex-determining factor sry",
										"sex-determining region y", "sry gene", "sry protein", "sex region of y chromosome");
		$synonym{$_} = "c-mos" foreach ("c-mos", "oocyte maturation factor mos", "mos", "mos protein");
		$synonym{$_} = "cryaa" foreach ("cryaa", "crya1", "alpha a-crystallin","alpha a-crystallin", "alphaa-crystallin",
										"alpha-a crystallin chain", "alpha-a-crystallin", "crystallin, alpha a", "alphaa-crystallin (crya1)",
										"alpha a crystallin");
		$synonym{$_} = "cryab" foreach ("cryab", "crya2", "alpha b-crystallin","alpha b-crystallin", "alphab-crystallin",
										"alpha-b crystallin chain", "alpha-b-crystallin", "crystallin, alpha b", "alphab-crystallin (crya2)",
										"alpha b crystallin");
		$synonym{$_} = "t-cell_receptor_beta" foreach ("t-cell receptor beta", "t cell receptor beta", "t-cell receptor beta chain",
										"t cell receptor beta chain", "t-cell receptor beta chain variable region",
										"t cell receptor beta chain variable region", "t-cell receptor beta chain variable",
										"t cell receptor beta chain variable", "t-cell receptor variable beta chain",
										"t cell receptor variable beta chain", "t-cell receptor beta-chain",
										"t cell receptor beta-chain", "t-cell receptor beta chain vj region",
										"t cell receptor beta chain vj region", "t-cell receptor beta chain v-d-j region",
										"t cell receptor beta chain v-d-j region", "t-cell receptor beta chain variable segment",
										"t cell receptor beta chain variable segment", "t-cell receptor beta chain constant region",
										"t cell receptor beta chain constant region", "t-cell receptor v beta gene",
										"t cell receptor v beta gene", "t-cell receptor-beta chain", "t cell receptor-beta chain",
										"t cell receptor bata chain", "t-cell receptor beta-chain", "t cell receptor beta-chain",
										"t-cell receptor beta chain (v-d-j-c)", "t cell receptor beta chain (v-d-j-c)",
										"t-cell receptor beta-chain v region", "t cell receptor beta-chain v region",
										"t-cell receptor v region beta-chain", "t cell receptor v region beta-chain",
										"t-cell receptor beta chain vdj region", "t cell receptor beta chain vdj region",
										"t-cell receptor beta chain v-region", "t cell receptor beta chain v-region",
										"t-cell receptor v-beta", "t cell receptor v-beta", "t-cell_receptor_beta");
		$synonym{$_} = "t-cell_receptor_alpha" foreach ("t-cell receptor alpha", "t cell receptor alpha", "t-cell receptor alpha chain",
										"t cell receptor alpha chain", "t-cell receptor alpha chain variable region",
										"t cell receptor alpha chain variable region", "t-cell receptor alpha chain variable",
										"t cell receptor alpha chain variable", "t-cell receptor variable alpha chain",
										"t cell receptor variable alpha chain", "t-cell receptor alpha-chain",
										"t cell receptor alpha-chain", "t-cell receptor alpha chain vj region",
										"t cell receptor alpha chain vj region", "t-cell receptor alpha chain v-d-j region",
										"t cell receptor alpha chain v-d-j region", "t-cell receptor alpha chain variable segment",
										"t cell receptor alpha chain variable segment", "t-cell receptor alpha chain constant region",
										"t cell receptor alpha chain constant region", "t-cell receptor v alpha gene",
										"t cell receptor v alpha gene", "t-cell receptor-alpha chain", "t cell receptor-alpha chain",
										"t cell receptor bata chain", "t-cell receptor alpha-chain", "t cell receptor alpha-chain",
										"t-cell receptor alpha chain (v-d-j-c)", "t cell receptor alpha chain (v-d-j-c)",
										"t-cell receptor alpha-chain v region", "t cell receptor alpha-chain v region",
										"t-cell receptor v region alpha-chain", "t cell receptor v region alpha-chain",
										"t-cell receptor alpha chain vdj region", "t cell receptor alpha chain vdj region",
										"t-cell receptor alpha chain v-region", "t cell receptor alpha chain v-region",
										"t-cell receptor v-alpha", "t cell receptor v-alpha", "t-cell_receptor_alpha");

	# tRNAs
		my @aminoAcids = qw(gly ala val leu ile met phe trp pro ser thr cys tyr asn gln asp glu lys arg his);
		my %AAcodes = ('glycine' => 'gly', 'alanine' => 'ala', 'valine' => 'val', 'leucine' => 'leu', 'isoleucine' => 'ile',
							'methionine' => 'met', 'phenylalanine' => 'phe', 'tryptophan' => 'trp', 'proline' => 'pro',
							'serine' => 'ser', 'threonine' => 'thr', 'cysteine' => 'cys', 'tyrosine' => 'tyr', 'asparagine' => 'asn',
							'glutamine' => 'gln', 'aspartic acid' => 'asp', 'glutamic acid' => 'glu', 'lysine' => 'lys',
							'arginine' => 'arg', 'histidine' => 'his');
		my %fullName = ('gly' => 'glycine', 'ala' => 'alanine', 'val' => 'valine', 'leu' => 'leucine', 'ile' => 'isoleucine',
							'met' => 'methionine', 'phe' => 'phenylalanine', 'trp' => 'tryptophan', 'pro' => 'proline',
							'ser' => 'serine', 'thr' => 'threonine', 'cys' => 'cysteine', 'tyr' => 'tyrosine', 'asn' => 'asparagine',
							'gln' => 'glutamine', 'asp' => 'aspartic acid', 'glu' => 'glutamic acid', 'lys' => 'lysine',
							'arg' => 'arginine', 'his' => 'histidine');

		foreach my $aa (@aminoAcids)
			{
			# Abbreviations
				$synonym{"trna $aa"} = "trna-$aa";
				$synonym{"transfer rna $aa"} = "trna-$aa";
				$synonym{"transfer rna-$aa"} = "trna-$aa";
				$synonym{"trna($aa)"} = "trna-$aa";
				$synonym{"trna $aa gene"} = "trna-$aa";
				$synonym{"trna-$aa gene"} = "trna-$aa";
				$synonym{"$aa trna"} = "trna-$aa";
				
				my $dashedName = $aa."-trna";
					$synonym{"$dashedName"} = "trna-$aa";
			
			# And again for full names
				$synonym{"trna $fullName{$aa}"} = "trna-$aa";
				$synonym{"transfer rna $fullName{$aa}"} = "trna-$aa";
				$synonym{"transfer rna-$fullName{$aa}"} = "trna-$aa";
				$synonym{"trna-$fullName{$aa}"} = "trna-$aa";
				$synonym{"trna($fullName{$aa})"} = "trna-$aa";
				$synonym{"trna $fullName{$aa} gene"} = "trna-$aa";
				$synonym{"trna-$fullName{$aa} gene"} = "trna-$aa";
				$synonym{"$fullName{$aa} trna"} = "trna-$aa";
				
				$dashedName = $fullName{$aa}."-trna";
				$synonym{"$dashedName"} = "trna-$aa";
			}
	}

sub geneClean
	{
	my $dirtyGene = shift;

	# Initial synonymizing
		$dirtyGene = lc($dirtyGene);	# Only work with lower-case for synonymizing and file names
		$dirtyGene = $synonym{$dirtyGene} if (defined $synonym{$dirtyGene});

	# Clean end of gene name
		$dirtyGene =~ s/\"$//g;
		$dirtyGene =~ s/\s+$//g;
		$dirtyGene =~ s/^\s+//g;
		
	# Remove punctuation that will confound file saving
		$dirtyGene =~ s/^'//g;
		$dirtyGene =~ s/'$//;
		$dirtyGene =~ s/:/ /g;
		$dirtyGene =~ s/, / /g;
		$dirtyGene =~ s/,//g;
		$dirtyGene =~ s/\*/-/g;
		$dirtyGene =~ s/\$/-/g;
		$dirtyGene =~ s/\#/-/g;
		$dirtyGene =~ s/\&/and/g;
		$dirtyGene =~ s/\//-/g;
		$dirtyGene =~ s/\\//g;
		$dirtyGene =~ s/\|/-/g;
		$dirtyGene =~ s/;//g;
		$dirtyGene =~ s/\<//g;
		$dirtyGene =~ s/\>//g;
		$dirtyGene =~ s/-+/-/g;
		$dirtyGene =~ s/\s+-/-/g;
		$dirtyGene =~ s/-\s+/-/g;
		$dirtyGene =~ s/^-+//;
		$dirtyGene =~ s/`/'/;
		
	# Collapse multiple whitespace
		$dirtyGene =~ s/\s+/_/g;
	
	# Clean up some tRNA variants (easier than specifying explicit synonyms for each tRNA)
		if ($dirtyGene =~ /tRNA/i)
			{
			$dirtyGene =~ s/_\d+$//;
			$dirtyGene =~ s/ \d+$//;
			}
	
	# Final synonymizing
		$dirtyGene = $synonym{$dirtyGene} if (defined $synonym{$dirtyGene});	# Recheck for synonym
		$dirtyGene = lc($dirtyGene);	# Ensure that gene names are in lower case for further synonymizing and file saving (final safety)
	
	return $dirtyGene;
	}

sub geneCount
	{
	# Set local parameters
		my $gbFile = shift;

		my (%iGene, @geneList, $speciesNum, %quickSpeciesGenePresent, %speciesCounter);
		my $wordBlock = 1;
		my $modelBlock = 0;
			my %modelSpecies;
		my $stripCount = 0;

	print "\nCounting gene occurrences in file $gbFile to establish genes that do not meet user-defined threshold(s)\n";

	# Search input file
		setLineBreak($gbFile);
		open (GENE, "<$gbFile") or die "Cannot open GenBank output file, $gbFile.\n";
			my $accCount = 0;
			my $modelFlag = 0;
			my $speciesFlag = 0;
			my $nameFlag = 0;
			my $countZero = time;
			my ($organism, $gene);
			while (<GENE>)
				{
				next unless ($_ =~ /LOCUS/ or $_ =~ /^\s*(\/gene)|(\/product)=\"/ or $_ =~ /^\s*\/organism=\"/ or $nameFlag == 1);
				chomp;
				my $gbLine = $_;
		
				if ($gbLine =~ /LOCUS/)
					{
					undef %iGene;
					undef $organism;
					$modelFlag = $speciesFlag = $nameFlag = 0;
					$accCount++;
					print "\t$accCount sequences read in\n" if ($accCount == int($accCount/10000)*10000 and $verbose);
					next;
					}
	
				# Get organism name
					if ($gbLine =~ /^\s*\/organism=\"(.+)\"/)
						{
						$organism = $1;
							$organism =~ s/\s/_/g;
						$modelFlag = 1 if (defined $modelSpecies{$organism});
						$speciesFlag = 1 if ($organism =~ /sp\./ or $organism =~ /cf\./);
						unless ($speciesFlag or ($modelFlag and $modelBlock))
							{
							$speciesCounter{$organism}++;
							$speciesNum++ if ($speciesCounter{$organism} == 1);
							}
						}
					next if ($modelFlag and $modelBlock);	# Entry pertains to model organism; skip parsing rest of entry
						
				# Get gene name
					if ($gbLine =~ /^\s*(\/gene)|(\/product)=\"(.+)/)	# Get start of gene name
						{
						$gene = $1 if ($gbLine =~ /=\"(.+)/);
						
						# Check whether name wraps onto a new line
							if (substr($gbLine, -1) ne "\"")
								{
								$nameFlag = 1;
								next;
								}
							else
								{
								$nameFlag = 2;
								$gene =~ s/\"$// if ($gene =~ /\"$/);
								}
						}
	
					if ($nameFlag == 1)	# Get continuation / end of gene name
						{
						$gene .= $gbLine;
						$nameFlag = 2 if ($gbLine =~ /\"$/)	# Gene name is complete
						}
	
					if ($gene and $nameFlag == 2)	# Process complete gene
						{
						next if (not defined $organism);
						unless ($wordBlock and ($gene =~ /hypothetical/i or $gene =~ /open reading frame/i or $gene =~ /similar/i or $gene =~ /homolog/i or $gene =~ /putative/i or $gene =~ /unknown/i or $gene =~ /unnamed/i or $gene =~ /\d+rik/i or $gene =~ /possible/i or $gene =~ /pseudo/i))
							{
							$gene = geneClean($gene);	# Clean up and synonymize gene name
		
							# Increment counters
								$iGene{$gene}++;	# Sequence counter
									$quickSequenceCount{$gene}++ if ($iGene{$gene} == 1);
	
								$quickSpeciesGenePresent{$gene}{$organism}++;	# Species-gene counter
									$quickSpeciesCount{$gene}++ if ($quickSpeciesGenePresent{$gene}{$organism} == 1);
							}
						$nameFlag = 0;
						undef $gene;
						}
				}
		close GENE;
		print "\n" if ($verbose);
		my $countTime = time - $countZero;
	
	# Print results and block genes as appropriate
		open (OUT, ">$countFile") or die "Cannot open results file for quick gene count, $countFile.\n";
			if ($speciesThreshold)
				{
				@geneList = sort {$quickSpeciesCount{$b} <=> $quickSpeciesCount{$a} } keys %quickSpeciesCount;
				}
			else
				{
				@geneList = sort {$quickSequenceCount{$b} <=> $quickSequenceCount{$a} } keys %quickSequenceCount;
				}

			printf OUT "Searching $accCount sequences took $countTime seconds with %s distinct genes found\n", scalar(@geneList);
				print OUT "\tModel organisms were excluded from search\n" if ($modelBlock);
				print OUT "\tCommon words were excluded from search\n" if ($wordBlock);
				print OUT "\nGene\tNo. of sequences\tNo. of species\n\n";
			
			print "\tGene counts:\n" if ($debug);
			foreach my $entry (@geneList)
				{
				print OUT "$entry\t$quickSequenceCount{$entry}\t$quickSpeciesCount{$entry}\n";
				print "\t$entry\t$quickSequenceCount{$entry}\t$quickSpeciesCount{$entry}\n" if ($debug);
				if ($quickSequenceCount{$entry} < $seqThreshold or $quickSpeciesCount{$entry} < $speciesThreshold)
					{
					$geneStatus{$entry} = "rejected";
					push @rejectList, $entry;
					}
				else
					{
					$stripCount++;
					}
				}
		close OUT;
		
		printf "\tSearching $accCount sequences took $countTime seconds with %s distinct genes found; results written to $countFile\n", scalar(@geneList);
		if (not $stripCount)
			{
			print "No genes were present more than the threshold value(s)\n";
			exit(0);
			}
	}

# Version history
#
#	v2.0 (July 25, 2005)
#		- improved and more complete parsing of GenBank files (esp. of gene names or CDs covering multiple lines,
#		  join statements and other multi-segment entries, detection of pseudogenes)
#		- added ability to mine input file for all gene content
#		- added ability to strip only those genes present for 1) a minimum number of sequences and/or species and
#		  2) species with valid names only
#		- added (incomplete) gene synonymy list (initially compiled by Robin Beck)
#		- minor bug fixes
#
#	v1.0.1a (April 29, 2004)
#		- added automatic detection of line breaks
#		- minor bug fixes
#
#	v1.0 (April 30, 2002)
#		- initial release
