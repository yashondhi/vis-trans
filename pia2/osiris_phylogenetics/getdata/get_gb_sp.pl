#!/usr/bin/env perl
use strict;
no warnings; #genbank produces annoying warning if no sequence is found

#use FindBin;
#use lib "$FindBin::Bin/lib";
use Bio::DB::GenBank;
use Bio::SeqIO;
use Bio::Root::Exception;
use Error qw(:try);


my $datafile = $ARGV[0];
my $datatype = $ARGV[1];
my $outtype = $ARGV[2];
my $outfile = $ARGV[3];
my $nodata = $ARGV[4];

my $accessions;
my @accnums;

	open (FILE,"<$datafile") or die "Cannot open file containing accession numbers\n";
	open (OUT,">$outfile") or die "Cannot open outfile\n";
	close OUT; #This overwrites old file if it exists
	open (ND,">$nodata") or die "Cannot open file\n";
	my $fh = Bio::SeqIO->newFh(-format=>$outtype, -file=>">>$outfile");


	while (<FILE>)
	{
	        chomp;
	        next unless ($_);
		push(@accnums, $_);
	}
	close FILE;

	my $countnames = 0;
	foreach (@accnums){
		#Should check input for one word per line and throw error if not, which is not done
	
	        $accessions = $_;
		chomp;
		if($accessions eq ""){
			die "Put spaces between accession numbers. No Empty Lines allowed.\n";
		}
	        my $qry_string .= $accessions."[organism]"." ";
	
#	        my $GBseq;
	        my $gb = new Bio::DB::GenBank;
	        my $query = Bio::DB::Query::GenBank->new
	                (-query   =>$qry_string,
	                 -db      =>$datatype);

		my $seqio;

		if (eval {$gb->get_Stream_by_query($query)}){
			$seqio = $gb->get_Stream_by_query($query);
			while( my $GBseq = $seqio->next_seq ) {
				my $sequence = $GBseq;   # read a sequence object
				print $fh $sequence; # write a sequence object
			}
		}else{
			print ND "$accessions\n";
		}
	}
exit;


