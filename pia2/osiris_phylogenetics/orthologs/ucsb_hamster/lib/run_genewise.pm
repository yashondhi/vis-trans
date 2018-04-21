package run_genewise;
use strict;
$ENV{'WISECONFIGDIR'} = "/home/osiris/galaxy-dist/tools/osiris/orthologs/ucsb_hamster/lib/wisecfg";
# this module runs genewise on a DNA sequence and a protein sequence
# and then allows to parse this result.
# the constructor creates an object containing a reference to an array
# containing the file content
1;
sub new {
    my $self_tmp = [];
    my $self;
    my ($class, $dna, $prot, $path) = @_;
    if (!defined $path) {
	$path = '/tmp';
    }
$dna =~ s/R/N/g;	#Added by THO -- genewise crashed with 'R' in dna sequence
$dna =~ s/S/N/g;	#Added by THO -- genewise crashed with 'R' in dna sequence
$dna =~ s/W/N/g;	#Added by THO -- genewise crashed with 'R' in dna sequence
$dna =~ s/D/N/g;	#Added by THO -- genewise crashed with 'R' in dna sequence
$dna =~ s/K/N/g;	#Added by THO -- genewise crashed with 'R' in dna sequence
$dna =~ s/Y/N/g;	#Added by THO -- genewise crashed with 'R' in dna sequence
$dna =~ s/B/N/g;	#Added by THO -- genewise crashed with 'R' in dna sequence
$dna =~ s/V/N/g;	#Added by THO -- genewise crashed with 'R' in dna sequence
$dna =~ s/M/N/g;	#Added by THO -- genewise crashed with 'R' in dna sequence

    # the file names
    my $protname = 'protein';
    my $dnaname = 'dna';
	
    ## print the two sequences to default path /tmp/
    open (DNA, ">$path/dna.fa") or die "could not open $path/dna.fa for writing\n";
    print DNA ">$dnaname\n$dna";
    close DNA;
    open (PROTEIN, ">$path/prot.fa") or die "could not open $path/prot.fa for writing\n";
    print PROTEIN ">$protname\n$prot";
    close PROTEIN;

    ## run genewise on the two sequences
  `echo \$WISECONFIGDIR`;
    
#    $self_tmp = [`.\/genewise -trans -cdna -pep -sum $path/prot.fa $path/dna.fa`];
#THO--For Galaxy run Genewise in the path
    $self_tmp = [`genewise -trans -cdna -pep -sum $path/prot.fa $path/dna.fa`];
    for (my $i = 0; $i < @$self_tmp; $i++) {
	$self_tmp->[$i] =~ s/\s{1,}$//;
    }
    $self->{gw} = $self_tmp;
    $self->{nt_seq} = $dna;
    $self->{prot_seq} = $prot;
    $self->{protname} = $protname;
    $self->{dnaname} = $dnaname;
    $self->{gw_count} = @$self_tmp;
    $self->{get_indel} = 1; ## per default the indel-part is recovererd, rather than masked by 'N'. See code for details
    $self->{indels} = _GetIndels($self_tmp);
    bless ($self, $class);
    return $self;}
#################
## sub score extract the score for the alignment
sub score {
    my $self = shift;
    my $score;
    for (my $i = 0; $i < $self->{gw_count}; $i ++) {
	if ($self->{gw}->[$i] =~ /^(\d{1,}\.{0,1}\d{0,}).*/) {
	    $score = $1;
	    last;
	}
    }
    return ($score);
}
##################
sub protein {
    my $self = shift;
    my $gw = $self->{gw};
    my $prot = '';
    for (my $i = 0; $i < @$gw; $i++) {
      if ($gw->[$i] =~ />.*\.pep/) { #the protein seq starts
	my $count = 1;
	while ($gw->[$i+$count] ne '//') {
	  my $protpart = $gw->[$i+$count];
	  chomp $protpart;
	  $prot .= $protpart;
	  $count ++;
	}
      }
      elsif (length $prot > 0) {
	last;
      }
    }
    return($prot);
 }
##################
sub translation {
    my $self = shift;
    my $finish = 0;
    my $translated_seq = '';
    my @transtmp;

    ## step 1: extract the relevant info from the genewise output
		
    for (my $i = 0; $i < $self->{gw_count}; $i++) {
      if ($self->{gw}->[$i] =~ />.*.tr/) {# a translated bit starts
	while ($self->{gw}->[$i] !~ '//') {
	  push @transtmp, $self->{gw}->[$i];
	  $i++;
	}
	last; # end the for loop since nothing left to be done
      }
    }
    
    ## step two: get the sequences
    my $count = -1;
    my $trans;
    for (my $i = 0; $i < @transtmp; $i++) {
      if ($transtmp[$i] =~ />/) {
	$count++;
	$trans->[$count]->{seq} = ''; # initialize
	if ($transtmp[$i] =~ /.*\[([0-9]{1,}):([0-9]{1,})\].*/) {
	  $trans->[$count]->{start} = $1;
	  $trans->[$count]->{end} = $2;
	  }
      }
      else {
	$trans->[$count]->{seq} .= $transtmp[$i];
      }
    }

    ## step 3: connect the fragments
    if (@$trans == 1) {
      $translated_seq = $trans->[0]->{seq};
    }
    else {
      for (my $i = 0; $i < @$trans; $i++) {
	$translated_seq .= $trans->[$i]->{seq};
	if ($i < (@$trans - 1)) {
	  my $missing = $trans->[$i+1]->{start} - $trans->[$i]->{end} -1;
	  $translated_seq .= 'X';
	}
      }
    }
    return($translated_seq);
  }

##################
sub codons {
    my $self = shift;
    my $finish = 0;
    my $codon_seq = '';
    my @transtmp;

    ## step 1: extract the relevant info from the genewise output
    for (my $i = 0; $i < $self->{gw_count}; $i++) {
      if ($self->{gw}->[$i] =~ />.*sp$/) {# the codons set starts
	while ($self->{gw}->[$i] !~ '//') {
	  push @transtmp, $self->{gw}->[$i];
	  $i++;
	}
	last; # end the for loop since nothing left to be done
      }
    }
    
    ## step two: get the sequences
    my $count = -1;
    my $trans;
    for (my $i = 0; $i < @transtmp; $i++) {
      if ($transtmp[$i] =~ />/) {
	$count++;
	$trans->[$count]->{seq} = ''; # initialize
	if ($transtmp[$i] =~ /.*\[([0-9]{1,}):([0-9]{1,})\].*/) {
	  $trans->[$count]->{start} = $1;
	  $trans->[$count]->{end} = $2;
	  }
      }
      else {
	$transtmp[$i] =~ tr/a-z/A-Z/;
	$trans->[$count]->{seq} .= $transtmp[$i];
      }
    }

    ## step 3: connect the fragments
    if ( @$trans == 1) {
      $codon_seq = $trans->[0]->{seq};
    }
    else {
      for (my $i = 0; $i < @$trans; $i++) {
	$codon_seq .= $trans->[$i]->{seq};
	if ($i < (@$trans - 1)) {
	  my $indel = '';
	  my $missing = $trans->[$i+1]->{start} - $trans->[$i]->{end} -1;

	  ## now decide whether the nts that did not got translated are masked by
	  ## 'N' or whether they will be represented as lower case letters
	  if ($self->{get_indel}) {
	    $indel = substr($self->{nt_seq}, $trans->[$i]->{end}, $missing);
	    $indel =~ tr/A-Z/a-z/;
	  }
	  else {
	    $indel = 'N' x $missing;
	  }
	  ## now append gap characters until the frame is recovered. Not that the gap
	  ## characters are added to the end of the indel-part. Thus, the codons are
	  ## not considered.
	  while (length($indel)%3 != 0) {
	    $indel .= '-';
	  }

	  $codon_seq .= $indel;
	}
      }
    }
    return ($codon_seq);
  }
###########################
sub protein_borders {
  my $self = shift;
  my $gw = $self->{gw};
  for (my $i = 0; $i < @$gw; $i++) {
    if ($gw->[$i] =~ /Bits.*introns$/) {
      my ($start, $end) = $gw->[$i+1] =~ /.*$self->{protname}\s{1,}([0-9]{1,})\s{1,}([0-9]{1,}).*/;
      return($start, $end);
    }
    else {
      die "no protein-start and end could not be determnined. Check genewise command\n";
    }
  }
}
##########################
sub cdna_borders {
  my $self = shift;
  my $gw = $self->{gw};
  for (my $i = 0; $i < @$gw; $i++) {
    if ($gw->[$i] =~ /Bits.*introns$/) {
      my ($start, $end) = $gw->[$i+1] =~ /.*$self->{dnaname}\s{1,}([0-9]{1,})\s{1,}([0-9]{1,}).*/;
      return($start, $end);
    }
    else {
      die "no cdna-start and end could not be determnined. Check genewise command\n";
    }
  }
}
##########################
sub _GetIndels {
  my $gw = shift;
  my $indel;
  for (my $i = 0; $i < @$gw; $i++) {
    if ($gw->[$i] =~ /Bits/) {
      $indel = $gw->[$i+1] =~ /.*([0-9]{1,})/;
      return($indel);
    }
  }
}
