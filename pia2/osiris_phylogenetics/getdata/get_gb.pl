#!/usr/bin/env perl
use strict;

#use FindBin;
#use lib "$FindBin::Bin/lib";
use Bio::DB::GenBank;
use Bio::SeqIO;


my $datafile = $ARGV[0];
my $datatype = $ARGV[1];
my $outtype = $ARGV[2];
my $outfile = $ARGV[3];
my $manual = $ARGV[4];
my $mannames = $ARGV[5];
my $genenames = $ARGV[6];


my $accessions;
my @accnums;
my @newnames;
my $manbin=0;
my @genenames;
my $genebin=0;

unless($mannames eq ''){
	@newnames = split(/ /,$mannames);
	$manbin=1;
}

unless($genenames eq ''){
	@genenames = split(/ /,$genenames);
	$genebin=1;
}

if($datafile eq 'None'){
	@accnums = split(/ /,$manual);
#	if(@accnums != @newnames && $manbin ==1 ){
#		die "Must have the same number of Custom Names as Accession Numbers\n";
#	}
}else{
	open (FILE,"<$datafile") or die "Cannot open file containing accession numbers\n";

	while (<FILE>)
	{
	        chomp;
	        next unless ($_);
		push(@accnums, $_);
	}
}
	my $countnames = 0;
	foreach (@accnums){
		#Should check input for one word per line and throw error if not, which is not done
	
	        $accessions = $_;
		chomp;
		if($accessions eq ""){
			die "Put spaces between accession numbers\n";
		}
	        my $qry_string .= $accessions."[accession]"." ";
	
	        my $GBseq;
	        my $gb = new Bio::DB::GenBank;
	        my $query = Bio::DB::Query::GenBank->new
	                (-query   =>$qry_string,
	                 -db      =>$datatype);
	
	        my $count;
	        my $species;
		my $seqio;
		if($outtype eq "phytab"){ #print phytab format, do not use bioperl as below.
			my $strand = 0;
			open(OUTFILE, ">>$outfile");
			if( defined ($seqio = $gb->get_Stream_by_query($query)) ){
			#        	my $seqio = $gb->get_Stream_by_query($query);
	        		while( defined ($GBseq = $seqio->next_seq )) {
	        		        my $sequence = $GBseq;   # read a sequence object
					for my $feat_object ($sequence->get_SeqFeatures) {          
					   $strand = $feat_object->strand;          
					}
					if($manbin ==1){ #replace GenBank Names with Custom Names
						$sequence->id($newnames[$countnames]);
						$sequence->desc('');
						$species = $sequence->id;
						$countnames++;
					}else{
						$species = $sequence->species->binomial;
						$species =~ s/ /_/g ;
					}
					if(@genenames > 0){
						if(@genenames == 1){
							if($strand < 0){
								my $revseq = $sequence->revcom;
		        					print OUTFILE $species."\t".$genenames[0]."\t".$sequence->accession."_REVCOMP\t".$revseq->seq."\n";
							}else{
		        					print OUTFILE $species."\t".$genenames[0]."\t".$sequence->accession."\t".$sequence->seq."\n";
							}
						}else{
							if($strand < 0){
								my $revseq = $sequence->revcom;
		        					print OUTFILE $species."\t".$genenames[0]."\t".$sequence->accession."_REVCOMP\t".$revseq->seq."\n";
							}else{
		        					print OUTFILE $species."\t".$genenames[0]."\t".$sequence->accession."\t".$sequence->seq."\n";
							}
					        }
					}else{
						if($strand < 0){
							my $revseq = $sequence->revcom;
		        				print OUTFILE $species."\tNone\t".$sequence->accession."_REVCOMP\t".$revseq->seq."\n";
						}else{
		        				print OUTFILE $species."\tNone\t".$sequence->accession."\t".$sequence->seq."\n";
						}
					}
	        		}
			}else{
				print "Did not find $accessions\n";
			}
		}else{
	        	my $fh = Bio::SeqIO->newFh(-format=>$outtype, -file=>">>$outfile");

			if( defined ($seqio = $gb->get_Stream_by_query($query)) ){
			#        	my $seqio = $gb->get_Stream_by_query($query);
	        		while( defined ($GBseq = $seqio->next_seq )) {
	        		        my $sequence = $GBseq;   # read a sequence object
					if($manbin ==1){ #replace GenBank Names with Custom Names
						$sequence->id($newnames[$countnames]);
						$sequence->desc('');
						$countnames++;
					}
	        		        print $fh $sequence; # write a sequence object
	        		}
			}else{
				print "Did not find $accessions\n";
			}
		}
	}
exit;
