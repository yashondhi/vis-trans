#!/usr/bin/perl -w
use strict;

use FindBin;
use lib "$FindBin::Bin/lib";
use Getopt::Long;
use Bio::SearchIO;
use Bio::Search::Hit::BlastHit;
use run_genewise;

# PROGRAMNAME: hamstrsearch_local.pl

# Copyright (C) 2009 INGO EBERSBERGER, ingo.ebersberger@univie.ac.at
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published
# by the Free Software Foundation; either version 3 of the License
# or any later version.

# This program is distributed in the hope that it will be useful
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program; If not, see http://www.gnu.org/licenses
 
# PROGRAM DESCRIPTION: This is the relevant version!

# DATE: Wed Dec 19 10:41:09 CEST 2007


######################## start main #############################
my $version = "\nhamstrsearch_local-v1.0.pl\nGalaxy implementation by Todd H. Oakley and Roger Ngo, UCSB\n";
print "$version\n";
### EDIT THE FOLLOWING LINES TO CUSTOMIZE YOUR SCRIPT
## note: all the variables can also be set via the command line
my $pid = $$;
my $prog = 'hmmsearch'; #program for the hmm search
my $alignmentprog = 'clustalw';


#PATH SETTINGS
my $hmmpath = '.';
my $blastpath = '.'; #Uses Galaxy working directory
my $tmpdir = 'tmp_' . $pid;
my $eval = 1; # eval cutoff for the hmm search
my $logfile = "hamstrsearch_" . $pid . '.log';
my $hmm_dir = 'hmm_dir';
my $fa_dir  = '';
##############################
my $help;
my $seq2store_file='';
my $cds2store_file='';
my $hmm;
my @hmms;
my $fa;
my $fafile;
my @seqs2store;
my @cds2store;
my $ep2eg;
my $estfile;
my $aln;
my $idfile;
my $taxon_check = 0;
my $hmmset;
my $hmmsearch_dir;
my $dbfile = ''; # the file hmmsearch is run against
my $dbfile_short;
my $taxon_file;
my $refspec_string;
my @refspec = qw();
my @primer_taxa;
my $refspec_name = '';
my $taxon_global;
my $fileobj;
my $fa_dir_neu = '';
my $gwrefprot;
my $seqtype;
my $align;
my $rep;
my $estflag;
my $proteinflag;
my $refseq;
my $refspec_final = '';
my $concat;
my $seqs2store_file;

#####determine the hostname#######
my $hostname = `hostname`;
chomp $hostname;
print "hostname is $hostname\n";

#################################
if (@ARGV==0) {
	$help = 1;
}
## help message
my $helpmessage = "
This program is freely distributed under a GPL. See -version for more info
Copyright (c) GRL limited: portions of the code are from separate copyrights

\nUSAGE: hamstrsearch_local.pl -sequence_file=<> -hmmset=<> -taxon=<>  -refspec=<> [-est|-protein] [-hmm=<>] [-representative] [-h]

OPTIONS:

-sequence_file: 
        path and name of the file containing the sequences hmmer is run against
	Per default, this file should be in the data directory.
-est:
        set this flag if you are searching in ESTs
-protein: 
        set this flag if you are searching in protein sequences
-hmmset: 
        specifies the name of the core-ortholog set.
        The program will look for the files in the default directory 'core-orthologs' unless you specify
	a different one.
-taxon:
        You need to specify a default taxon name from which your ESTs or protein sequences are derived.
-refspec: 
        sets the reference species. Note, it has to be a species that contributed sequences 
        to the hmms you are using. NO DEFAULT IS SET! For a list of possible reference
	taxa you can have a look at the speclist.txt file in the default core-ortholog sets
	that come with this distribution. Please use the 5 letter abreviations. If you choose
	to use core-orthologs were not every taxon is represented in all core-orthologs, you
	can provide a comma-separated list with the preferred refspec first. The lower-ranking 
        reference species will only be used if a certain gene is not present in the preferred 
        refspecies due to alternative paths in the transitive closure to define
	the core-orthologs.
        CURRENTLY NO CHECK IS IMPLEMENTED!
        NOTE: A BLAST-DB FOR THE REFERENCE SPECIES IS REQUIRED!
-eval_limit=<>
        This options allows to set the e-value cut-off for the HMM search.
        DEFAULT: 1
-hmm: 
        option to provide only a single hmm to be used for the search. 
        Note, this file has to end with .hmm 

### the following options should only be used when you chose to alter the default structure of the
### hamstrsearch_local directories. Currently, this has not been extensively tested.
-fasta_file: 
        path and name of the file containing the core-ortholog sequences
	you don't have to use this option when you 
-hmmpath: 
        sets the path to the hmm_dir. By default this is set to the current directory.
-blastpath: 
        sets the path where the blast-dbs are located. Per default this is ../blast_dir
        Note, the program expects the structure blastpath/refspec/refspec_prot.
        To overrule this structure you will have to edit the script.
        \n\n";
GetOptions ("h"        => \$help,
            "hmm=s"    => \$hmm,
            "est"    => \$estflag,
            "protein"=> \$proteinflag,
            "sequence_file=s" => \$dbfile,
            "fasta_file=s" => \$fafile,
            "hmmset=s" => \$hmmset,
            "hmmpath=s" => \$hmmpath,
            "taxon_file=s" => \$taxon_file,
            "taxon=s"  => \$taxon_global,
            "eval_limit=s" => \$eval,
            "refspec=s" => \$refspec_string,
            "estfile=s" => \$estfile,
            "representative" => \$rep,
            "blastpath=s" => \$blastpath,
	    "galaxyout=s" => \$seqs2store_file,
	    "2galaxyout=s" => \$cds2store_file);

if ($help) {
  print $helpmessage;
  exit;
}

## 1) check if all information is available to run HaMStR
my ($check, @log) = &checkInput();
if ($check == 0) {
  print join "\n", @log;
  print "$helpmessage";
  exit;
}
else {
  open (OUT, ">$logfile") or die "could not open logfile $logfile\n";
  print OUT join "\n", @log;
  close OUT;
}
### read in of the core-ortholog sequences
my $co_seqs = parseSeqfile("$fafile");

## 2) loop through the hmms
## process each hmm file separately
for (my $i = 0; $i < @hmms; $i++) {
  $fileobj = undef;
  my @seqs = qw();
  my @newseqs = qw();## var to contain the sequences to be added to the orthologous cluster
  my @newcds = qw();
  my $hmm = $hmms[$i];
  my $hmmout = $hmm;
  $hmmout =~ s/\.hmm/\.out/;
  ## 3) run the hmm search
  if (!(-e "$hmmsearch_dir/$hmmout")) {
    print "now running $prog using $hmm\n";
    !`$prog $hmm_dir/$hmm $dbfile >$hmmsearch_dir/$hmmout` or die "Problem running hmmsearch\n";
  }
  else {
    print "an hmmresult $hmmout already exists. Using this one!\n";
    print "NOTE: in Galaxy the hmm results are stored in the directory of the dataset\n";
  }
  
  ## 4) process the hmm search result
  my $hitcount = 0;
  ## 4a) loop through the individual results
  ## now the modified version for hmmer3 comes
  my ($query_name, @results) = parseHmmer3($hmmout, $hmmsearch_dir);
  if (! @results) {
    print "no hit found for $query_name\n";
    next;
  }
  chomp $query_name;
  print "Results for $query_name\n";
  for (my $k = 0; $k < @results; $k++) {
    my $hitname = $results[$k];
    print "$hitname\n";
    my $keep = 0;
    my $hitseq = '';
    $refseq = '';
    ## 4b) test for the reciprocity criterion fulfilled
    ($keep, $hitname, $hitseq, $refspec_final, $refseq)  = &check4reciprocity($query_name, $hitname, @refspec);
    if ($keep == 1) {
      ## blast search with the hmm hit identifies the core-ortholog sequence of the reference species
      ## check for the taxon from the taxon_file.txt.
      my $taxon = '';
      if ($taxon_check){
	if ($taxon_check == 1) {
	  $taxon = &getTaxon($hitname);
	}
	elsif ($taxon_check == 2) {
	  $taxon = $taxon_global;
	}
      }
      ## put the info about the hits into an object for later post-processing
      $fileobj->{$taxon}->{prot}->[$hitcount] = $hitseq;
      $fileobj->{$taxon}->{ids}->[$hitcount] = $hitname;
      $fileobj->{$taxon}->{refseq} = $refseq;
      $hitcount++;
    }
    else {
      print "match to different protein from $refspec_final\n";
    }
  }
  ## 5) do the rest only if at least one hit was obtained
  if (defined $fileobj) {
    ## 5a) if the hits are derived from ESTs, get the best ORF
    if ($estflag) {
      $fileobj =  &predictORF();
    }
    ## 5b) if the user has chosen to postprocess the results
    if ($rep) {
      &processHits($refseq, $concat);
    }
    ## 6) prepare the output
    my @taxa = keys(%$fileobj);
    for (my $i = 0; $i< @taxa; $i++) {
      if ($rep) {
#	push @newseqs, ">$query_name|$taxa[$i]|$fileobj->{$taxa[$i]}->{refid}";
#Rearrange Order for Galaxy - want species first for phylotable format convert pipes to tabs
	push @newseqs, ">$taxa[$i]\t$query_name\t$fileobj->{$taxa[$i]}->{refid}";
	push @newseqs, $fileobj->{$taxa[$i]}->{refprot};
	if ($estflag) {
#	  push @newcds, ">$query_name|$taxa[$i]|$fileobj->{$taxa[$i]}->{refid}";
#Rearrange Order for Galaxy - want species first for phylotable format
	  push @newcds, ">$taxa[$i]\t$query_name\t$fileobj->{$taxa[$i]}->{refid}";
	  push @newcds, $fileobj->{$taxa[$i]}->{refcds};
	}
      }
      else {
	my $idobj = $fileobj->{$taxa[$i]}->{ids};
	my $protobj = $fileobj->{$taxa[$i]}->{prot};
	my $cdsobj  = $fileobj->{$taxa[$i]}->{cds};
	for (my $j = 0; $j < @$idobj; $j++) {
#	  push @newseqs, ">$query_name|$taxa[$i]|$idobj->[$j]";
#Rearrange Order for Galaxy - want species first for phylotable format also tabs not pipe
	  push @newseqs, ">$taxa[$i]\t$query_name\t$idobj->[$j]";
	  push @newseqs, $protobj->[$j];
	  if ($estflag) {
#	    push @newcds, ">$query_name|$taxa[$i]|$idobj->[$j]";
#Rearrange Order for Galaxy - want species first for phylotable format
	    push @newcds, ">$taxa[$i]\t$query_name\t$idobj->[$j]";
	    push @newcds, $cdsobj->[$j];
	  }
	}
      }
      my $refs = $co_seqs->{$query_name};
      for (keys %$refs) {
	my $line = ">$query_name|$_|" . $refs->{$_}->{seqid} . "\n" . $refs->{$_}->{seq};
	push @seqs, $line;
      }
      chomp @seqs;
      print "\n";
      @seqs = (@seqs, @newseqs);
      open (OUT, ">$fa_dir_neu/$query_name.fa");
      print OUT join "\n", @seqs;
      print OUT "\n";
      close OUT;
      for (my $i = 0; $i < @newseqs; $i+= 2) {
#	my $line = $newseqs[$i] . "|" . $newseqs[$i+1];
#Galaxy uses tabs not pipes
	my $line = $newseqs[$i] . "\t" . $newseqs[$i+1];
	$line =~ s/>//;
	push @seqs2store, $line;
	if ($estflag) {
#Galaxy uses tabs not pipes
#	  my $cdsline = $newcds[$i] . "|" . $newcds[$i+1];
	  my $cdsline = $newcds[$i] . "\t" . $newcds[$i+1];
	  $cdsline =~ s/>//;
	  push @cds2store, $cdsline;
	}
      }
    }
  }
}
if (@seqs2store > 0) {
  open (OUT, ">$seqs2store_file") or die "failed to open output SEQS file\n";
  print OUT join "\n", @seqs2store;
  print OUT "\n";
  close OUT;
  if ($estflag) {
    open (OUT, ">$cds2store_file") or die "failed to open output CDS file\n";
    print OUT join "\n", @cds2store;
    print OUT "\n";
    close OUT;
  }
}
else {
  open (OUT, ">$seqs2store_file") or die "failed to open output SEQS file\n";
  print OUT "no hits found\n";
}
exit;
##################### start sub ###############
####### checkInput performs a number of checks whether sufficient information
### and all data are available to run HaMStR
sub checkInput {
  my @log;
  my $check = 1;
  $dbfile_short = $dbfile;
  $dbfile_short =~ s/\..*//;
  ## 1) check for filetype
  print "Checking for filetype:\t";
  if (!defined $estflag and !defined $proteinflag) {
    push @log, "please determine the sequence type. Choose between the options -EST or -protein";
    print "failed\n";
    $check = 0;
  }
  else {
    if ($estflag) {
      $estfile = $dbfile;
      $dbfile = "$dbfile.tc";
      push @log, "HaMStR will run on the ESTs in $estfile";
      push @log, "Translating ESTs";
      if (!(-e "$dbfile")) {
	print "translating $estfile, this may take a while\n";
  	`translate.pl -in=$estfile -out=$dbfile`;
	open (LOG, "$logfile") or die "could not open logfile $logfile\n";
	my @info = <LOG>;
	@log = (@log, @info);
	close LOG;
      }
      else {
	push @log, "Translated file already exists, using this one\n";
      }
      if (! -e "$dbfile") {
	push @log, "The translation of $estfile failed. Check the script translate.pl";
	print "failed\n";
	$check = 0;
      }
    }
    else {
      ## file type is protein
      print "succeeded\n";
    }
  }
  ## 2) Check for presence of blastall
  print "Checking for the blast program\t";
  if (!(`blastall`)) {
    push @log, "could not execute blastall. Please check if this program is installed and executable";
    print "failed\n";
    $check = 0;
  }
  else {
    push @log, "check for blastall succeeded";
    print "succeeded\n";
  }
  ## 3) Check for presence of hmmsearch
  print "Checking for hmmsearch\t";
  if (! `$prog -h`) {
    push @log, "could not execute $prog. Please check if this program is installed and executable";
    print "failed\n";
    $check = 0;
  }
  else {
      push @log, "check for $prog succeeded\n";
      print "succeeded\n";
  }
  ## 4) Check for reference taxon
  print "Checking for reference species and blast-dbs\t";
  if (!(defined $refspec_string)) {
      push @log, "Please provide a reference species for the reblast!";
      print "failed\n";
      $check = 0;
  }
  else {
    push @log, "Reference species for the re-blast: $refspec_string";
    @refspec = split /,/, $refspec_string;
    $refspec_name = $refspec[0];
    print "succeeded\n";
  }
  ## 5) Check for presence of the required blast dbs
  print "checking for blast-dbs:\t";
  for (my $i = 0; $i < @refspec; $i++) {
    my $blastpathtmp = "$blastpath/$refspec[$i]/$refspec[$i]" . "_prot";
    if (! (-e "$blastpathtmp.pin")) {
      push @log, "please edit the blastpath. Could not find $blastpathtmp";
      print "$blastpathtmp failed\n";
      $check = 0;
    }
    else {
      push @log, "check for $blastpathtmp succeeded";
      print "succeeded\n";
    }
  }
  ## 6) Check for presence of the directory structure
  print "checking for presence of the hmm files:\t";
  if (!(defined $hmmset)) {
    $hmmpath = '.';
    $hmmset = 'manual';
  }
  else {
    $hmmpath = "$hmmpath/$hmmset";
    $fafile = "$hmmpath/$hmmset" . '.fa';
  }
  $hmm_dir = "$hmmpath/$hmm_dir";
  $hmmsearch_dir = $dbfile_short . '_' . $hmmset;

#CHANGED FOR GALAXY DIRECTORY
 # $fa_dir_neu = 'fa_dir_' . $dbfile_short . '_' . $hmmset . '_' . $refspec_name;
  $fa_dir_neu = $dbfile_short . '_' . $hmmset . '_' . $refspec_name;
  ## 7) check for the presence of the hmm-files and the fasta-file
  if (!(-e "$hmm_dir")) {
    push @log, "Could not find $hmm_dir";
    print "failed\n";
    $check = 0;
  }
  else {
    if (defined $hmm) {
      if (! -e "$hmm_dir/$hmm") {
	push @log, "$hmm has been defined but could not be found in $hmm_dir/$hmm";
	$check = 0;
      }
      else {
	push @log, "$hmm has been found";
	if ($hmm =~ /\.hmm$/) {
	  @hmms = ($hmm);
	}
      }
    }
    else {
      push @log, "running HaMStR with all hmms in $hmm_dir";
      @hmms = `ls $hmm_dir`;
    }
    chomp @hmms;
    print "succeeded\n";
  }
  
  ## 8) Test for presence of the fasta file containing the sequences of the core-ortholog cluster
  print "checking for presence of the core-ortholog file:\t";
  if (defined $fafile) {
    if (! -e "$fafile") {
      push @log, "Could not find the file $fafile";
      print "failed\n";
      $check = 0;
    }
    else {
      push @log, "check for $fafile succeeded";
      print "succeeded\n";
    }
  }
  else {
    push @log, "Please provide path and name of fasta file containing the core-ortholog sequences";
    $check = 0;
    print "failed\n";
  }
  ## 9) Checks for the taxon_file
  print "testing whether the taxon has been determined:\t";  
  if (!(defined $taxon_file) or (!(-e "$taxon_file"))) {
    if (defined $taxon_global) {
      push @log, "using default taxon $taxon_global for all sequences";
      print "succeeded\n";
      $taxon_check = 2;
    }
    else {
      push @log, "No taxon_file found. Please provide a global taxon name using the option -taxon";
      print "failed\n";
      $check = 0;
    }
  }
  else {
    push @log, "using the file $taxon_file as taxon_file";
    print "succeeded\n";
    $taxon_check = 1;
  }
  ## 10) Set the file where the matched seqs are found
#CHANGED BY THO FOR GALAXY TO ALLOW DETERMINATION OF OUTPUT FILE Made INPUT Option
#  $seqs2store_file = 'hamstrsearch_' . $dbfile_short . '_' . $hmmset . '.out';
#  $cds2store_file = 'hamstrsearch_' . $dbfile_short . '_' . $hmmset . '_cds.out';

  ## 11) apply the evalue-cut-off to the hmmsearch program
  $prog = $prog . " -E $eval";
  push @log, "hmmsearch: $prog";

  ## 12) setting up the directories where the output files will be put into.
  if ($check == 1) {
    if (!(-e "$hmmsearch_dir")) {
      `mkdir $hmmsearch_dir`;
    }
    if (!(-e "$fa_dir_neu")) {
      `mkdir $fa_dir_neu`;
    }
    if (!(-e "$tmpdir")) {
      `mkdir $tmpdir`;
    }
  }
  return ($check, @log);
}
#################
## check4reciprocity is the second major part of the program. It checks
## whether the protein sequence that has been identified by the hmmsearch
## identifies in turn the protein from the reference taxon that was used to
## build the hmm.
sub check4reciprocity {
  my ($query_name, $hitname, @refspec) = @_;
  my $searchdb;
  ## get the sequence from the db_file
  my $hitseq = `grep -m 1 -A 1 ">$hitname" $dbfile | tail -n 1`;
  if (!defined $hitseq) {
    print "could not retrieve a sequence for $hitname. Skipping...\n";
    return(0, '', '', '');
  }
  else {
    ## now get the sequence used for building the hmm
    my @original;
    my $refspec_final;
    for (my $i = 0; $i < @refspec; $i++) {
      @original = `grep "^>$query_name|$refspec[$i]" $fafile |sed -e "s/.*$refspec[$i]\|//"`;
      chomp @original;
      
      if (@original > 0) {
	$refspec_final = $refspec[$i];
	$searchdb = "$blastpath/$refspec_final/$refspec_final" . "_prot";
	last;
      }
      else {
	print "original sequence not be found with grepping for ^>$query_name|$refspec[$i]. Proceeding with next refspec\n";
      }
    }
    if (@original == 0) {
      print "original sequence not be found\n";
      return (0, '', '', $refspec_final);
    }
    print "REFSPEC is $refspec_final\n";
    ## continue with the blast
    chomp $hitseq;
#    $hitname =~ s/\|.*//;
    ## now run the blast
    open (OUT, ">$tmpdir/$$.fa") or die "could not open out for writing\n";
    print OUT ">$hitname\n$hitseq";
    close OUT;
    !`blastall -p blastp -d $searchdb -v 10 -b 10 -i $tmpdir/$$.fa -o $tmpdir/$$.blast` or die "Problem running blast\n";
    ## now parse the best blast hit
    my @hits = &getBestBlasthit("$tmpdir/$$.blast");
    if (@hits > 0) {
      
      print "hmm-seq: ", join "\t", @original , "\n";
      ## now loop through the best hits with the same evalue and check whether
      ## among these I find the same seq as in $original
      for (my $i = 0; $i <@hits; $i++) {
	print "blast-hit: $hits[$i]";
	## now loop through all the refspec-sequences in the hmm file
	for (my $j = 0; $j < @original; $j++) {
	  if ($original[$j] eq $hits[$i]) {
	    print "\tHit\n";
	    my ($refseq) = `grep -A 1 "$query_name|$refspec_final|$original[$j]" $fafile |tail -n 1`;
	    return (1, $hitname, $hitseq, $refspec_final, $refseq);
	  }
	  else {
	    print "\nnot hitting $original[$j]\n";
	  }
	}
      }
      ### if we end up here, we didn't find a hit that matches to the original sequence
      ### in the top hits with the same eval
      return (0, '', '', $refspec_final);
    }
    else {
      print "no hit obtained\n";
      return(0, '', '', $refspec_final);
    }	
  }
}
#############
sub getBestBlasthit {
    my @hits;
    my ($file) = @_;
    my $searchio = Bio::SearchIO->new(-file        => $file,
				      -format      => 'blast',
				      -report_type => 'blastp') or die "parse failed";
    while( my $result = $searchio->next_result ){
	my $count = 0;
	my $sig;
	my $sig_old;
	while( my $hit = $result->next_hit){
	    ## now I enter all top hits having the same evalue into the result
	    $sig = $hit->score;
	    if (!defined $sig_old) {
		$sig_old = $sig;
	    }
	    if ($sig == $sig_old) {
		push @hits, $hit->accession;
	    }
	    else {
		last;
	    }
	}
    }
    return(@hits);
}
##################
sub getTaxon {
    my ($hitname) = @_;
#    my $q = "select name from taxon t, est_project e, est_info i, annotation_neu a where a.id = $hitname and a.contig_id = i.contig_id and i.project_id = e.project_id and e.taxon_id = t.taxon_id";
    if ($hitname =~ /\D/) {
	$hitname =~ s/_.*//;
    }
    my $taxon = `grep -m 1 "^$hitname," $taxon_file | sed -e 's/^.*,//'`;
    chomp $taxon;
    $taxon =~ s/^[0-9]+,//;
    $taxon =~ s/\s*$//;
    $taxon =~ s/\s/_/g;
    if ($taxon) {
	return ($taxon);
    }
    else {
	return();
    }
}
###############
sub processHits {
  my ($concat) = @_; 
  ## 1) align all hit sequences for a taxon against the reference species
  my @taxa = keys(%$fileobj);
  for (my $i = 0; $i < @taxa; $i++) {
    &orfRanking($taxa[$i]);
  }
}  

################
sub predictORF {
  my $fileobj_new;
#  my ($refseq) = @_;
  my @taxa = keys(%$fileobj);
  for (my $i = 0; $i < @taxa; $i++) {
    my $protobj = $fileobj->{$taxa[$i]}->{prot};
    my $idobj = $fileobj->{$taxa[$i]}->{ids};
    my $refseq = $fileobj->{$taxa[$i]}->{refseq};
    my @ids = @$idobj;
    for (my $j = 0; $j < @ids; $j++) {
      ## determine the reading frame
      my ($rf) = $ids[$j] =~ /.*_RF([0-9]+)/;
	print "rf is $rf\n";
      $ids[$j] =~ s/_RF.*//;
#      my $est = `grep -A 1 "$ids[$j]" $estfile |tail -n 1`;
################new grep command from version 8 to fix bug   
       my $est = `grep -A 1 ">$ids[$j]\\b" $estfile |tail -n 1`;
      if (! $est) {
	die "error in retrieval of est sequence for $ids[$j] in subroutine processHits\n";
      }
      ## the EST is translated in rev complement
      if ($rf > 3) {
	$est = revComp($est);
      }

      my $gw = run_genewise->new($est, $refseq, "$tmpdir");
      my $translation = $gw->translation;
      my $cds = $gw->codons;
      $translation =~ s/[-!]//g;
      $fileobj_new->{$taxa[$i]}->{ids}->[$j] = $ids[$j];
      $fileobj_new->{$taxa[$i]}->{prot}->[$j] = $translation;
      $fileobj_new->{$taxa[$i]}->{cds}->[$j] = $cds;
      $fileobj_new->{$taxa[$i]}->{refseq} = $refseq;
    }
  }
  return($fileobj_new);
}
############################
sub orfRanking {
  my ($spec) = @_;
  my $result;
  my $refprot;
  my $refcds;
  my @toalign;
  my $protobj = $fileobj->{$spec}->{prot};
  my $idobj = $fileobj->{$spec}->{ids};
  my $refcluster; ## variables to take the cluster and its id for later analysis
  my $refid;
  if (@$protobj == 1) {
    ## nothing to chose from
    $refprot = $protobj->[0];
    $refcds = $fileobj->{$spec}->{cds}->[0];
    my $length = length($refprot);
    $refid = $idobj->[0] . "-" . $length;
  }
  else {
    ## more than one cluster
    push @toalign, ">$refspec_final";
    push @toalign, $fileobj->{$spec}->{refseq};
    ## now walk through all the contigs
    for (my $i = 0; $i < @$protobj; $i++) {
      my @testseq = (">$idobj->[$i]", $protobj->[$i]);
      @testseq = (@testseq, @toalign);
      open (OUT, ">$tmpdir/$pid.ref.fa") or die "could not open file for writing refseqs\n";
      print OUT join "\n", @testseq;
      close OUT;
      ## run clustalw
      !(`$alignmentprog $tmpdir/$pid.ref.fa -output=fasta -outfile=$tmpdir/$pid.ref.aln 2>&1 >$tmpdir/$pid.ref.log`) or die "error running clustalw\n";
      ## get the alignment score
      $result->[$i]->{score} =  `grep "Alignment Score" $tmpdir/$pid.ref.log |sed -e 's/[^0-9]//g'`;
      if (!$result->[$i]->{score}) {
	die "error in determining alignment score\n";
      }
      chomp $result->[$i]->{score};
      ## get the aligned sequence
      open (ALN, "$tmpdir/$pid.ref.aln") or die "failed to open alignment file\n";
      my @aln = <ALN>;
      close ALN;
      my $aseq = extractSeq($idobj->[$i], @aln);
      ## remove the terminal gaps
      $aseq =~ s/-*$//;
      $result->[$i]->{aend} = length $aseq;
      my ($head) = $aseq =~ /^(-*).*/;
      ($result->[$i]->{astart}) = length($head)+1;
    }
    ### the results for all seqs has been gathered, now order them
    $result = sortRef($result);
    ($refprot, $refcds, $refid) = &determineRef($result,$spec);
  }
  $fileobj->{$spec}->{refprot} = $refprot;
  $fileobj->{$spec}->{refcds}  = $refcds;
  $fileobj->{$spec}->{refid}   = $refid;
  return();
}
###########################
sub sortRef {
    my $result = shift;
    my @sort;
    for (my $i = 0; $i < @$result; $i++) {
	push @sort, "$i,$result->[$i]->{astart},$result->[$i]->{aend},$result->[$i]->{score}";
    }
    open (OUT, ">$tmpdir/$pid.sort") or die "failed to write for sorting\n";
    print OUT join "\n", @sort;
    close OUT;
    `sort -n -t ',' -k 2 $tmpdir/$pid.sort >$tmpdir/$pid.sort.out`;
    @sort = `less $tmpdir/$pid.sort`;
    chomp @sort;
    $result = undef;
    for (my $i = 0; $i < @sort; $i++) {
	($result->[$i]->{id}, $result->[$i]->{start}, $result->[$i]->{end}, $result->[$i]->{score}) = split ',', $sort[$i];
    }
    return($result);
}
########################
sub determineRef {
  my ($result, $spec) = @_;
  my $lastend = 0;
  my $lastscore = 0;
  my $final;
  my $count = 0;
  my $id = '';
  for (my $i = 0; $i < @$result; $i++) {
    if ($result->[$i]->{start} < $lastend or $lastend == 0) {
      if ($result->[$i]->{score} > $lastscore) {
	$lastend = $result->[$i]->{end};
	$lastscore = $result->[$i]->{score};
	$id = $result->[$i]->{id};
      }
    }
    elsif ($result->[$i]->{start} > $lastend) {
      ## a new part of the alignment is covered. Fix the results obtained so far
      $final->[$count]->{id} = $id;
      $lastend = $result->[$i]->{end};
      $id = $result->[$i]->{id};
      $count++;
    }
  }
  $final->[$count]->{id} = $id;
  ## now concatenate the results
  my $refprot = '';
  my $refid = '';
  my $refcds = '';
  for (my $i = 0; $i < @$final; $i++) {
    my $seq = $fileobj->{$spec}->{prot}->[$final->[$i]->{id}];
    my $cdsseq = $fileobj->{$spec}->{cds}->[$final->[$i]->{id}];
    my $length = length($seq);
    $refid .= "$fileobj->{$spec}->{ids}->[$final->[$i]->{id}]-$length" . "PP";
    $refprot .= $seq;
    if ($estflag) {
      $refcds .= $cdsseq;
    }
  }
  $refid =~ s/PP$//;
  return($refprot, $refcds, $refid);
}
#############################
sub extractSeq {
  my ($id, @aln) = @_;
  my $seq = '';
  my $start = 0;
  for (my $i = 0; $i < @aln; $i++) {
    if ($aln[$i] =~ $id) {
      $start = 1;
    }
    elsif ($aln[$i] =~ />/ and $start == 1) {
      last;
    }
    elsif ($start == 1) {
      $seq .= $aln[$i];
    }
  }
  $seq =~ s/\s//g;
  return ($seq);
}
##############################
sub revComp {
    my ($seq) = @_;
    $seq =~ tr/AGCTYRKMWSagct/TCGARYMKWSTCGA/;
    $seq = reverse($seq);
    return($seq);
}
##############################
sub parseHmmer3 {
  my ($file, $path) = @_;
  if (!defined $path) {
    $path = '.';
  }
  open (IN, "$path/$file") or die "failed to open $file\n";
  my @data = <IN>;
  close IN;
  ### extract the hits
  my @hit;
  my $start = 0;
  my $stop = 0;
  my $i = 0;
  for (my $i = 0; $i < @data; $i++) {
    if (!($data[$i] =~ /\S/)) {
      next;
    }
    else {
      if ($data[$i] =~ /Scores for complete sequence/) {
	$start = 1;
	$i += 4;
      }
      elsif (($data[$i] =~ /inclusion threshold/) or ($data[$i] =~ /Domain and alignment/)) {
	last;
      }
      if ($start == 1 and $stop == 0) {
	$data[$i] =~ s/^\s+//;
	my @list = split /\s+/, $data[$i];
	push @hit, $list[8];
	$start = 0; #Added by THO
      }
    }
  }
  ### get the query_id
  my ($query) = grep /^Query:/, @data;
  $query =~ s/^Query:\s+//;
  $query =~ s/\s.*//;
  if (defined $hit[0]) {
    chomp @hit;
    return ($query, @hit);
  }
  else {
    return ($query);
  }
}
#####################
sub parseSeqfile {
  my $seqref;
  my $id;
  my $spec;
  my $seqid;
  my $seq;
  my $file = shift;
  open (IN, "$file") or die "failed to open $file\n";
  my @seqs = <IN>;
  close IN;
  chomp @seqs;
  for (my $i = 0; $i < @seqs; $i++) {
    if ($seqs[$i] =~ />/) {
	$seqs[$i] =~ s/>//;
      if (defined $id and defined $seq) {
	$seqref->{$id}->{$spec}->{seqid} = $seqid;
	$seqref->{$id}->{$spec}->{seq} = $seq;
	$seq = undef;
      }
      ($id, $spec, $seqid) = split (/\|/, $seqs[$i]);
    }
    else {
      $seq .= $seqs[$i];
    }
  }
  if (defined  $id and defined $seq) {
	$seqref->{$id}->{$spec}->{seqid} = $seqid;
	$seqref->{$id}->{$spec}->{seq} = $seq;
	$seq = undef;
      }
  return ($seqref);
}
