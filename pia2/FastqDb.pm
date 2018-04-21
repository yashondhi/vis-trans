=head1 NAME

Fastq - Simple object for Fastq sequence

=head1 SYNOPSIS

    my $seq=new Fastq( $hdr, $id, $base, $barcode, $pair, $seq, $qual, $convert_qual_to_sanger );
    $seq->qc;
    print $seq->output;

=head1 DESCRIPTION

Object for a single read sequence, with methods for basic manipulation and quality control.

=head1 METHODS

=over 5

=cut

package Fastq;

use constant {
    CHARACTERS_PER_LINE => 80, # for formatting Fastq output
    CLOSE_ENOUGH_TO_END => 6,   # hits of adapters this close to end will result in trimming to be done to end
    DEFAULT_MINLEN => 20,
    DEFAULT_MEANQ => 20,
    DEFAULT_WINSIZE => 5,
    DEFAULT_LOW_COMPLEXITY => 0.8,
    DEFAULT_MAXN => 3,
};

=item new $hdr $seq $qual $convert

Initialize new sequence object. If the quality scores use Illumina scaling, the $convert flag *must* be set as the object
assumes and requires sanger-scaling.  If initialized with Illumina-scaled Q-scores and the convert flag was not set,
do not use any methods which use Q-scores (e.g. qual_end_trim, output_qual).

=cut

sub new {
    my ($class,$hdr,$seq,$qual,$params)=@_;
    die("Incomplete Fastq record: hdr=$hdr, seq=$seq, qual=$qual\n") unless $hdr and $seq and $qual;
    # INIT
    my $this={
        seq=>$seq, # complete sequence without newlines
        qual=>$qual # complete quality string without newlines
    };
    bless $this,$class;
    $this->_parse_header($hdr,$params->{format}); # populates id, base, pair, barcode values
    $this->_trim_roche_mid($params->{roche_mid_len}) if exists($params->{roche_mid_len});
    $this->convert_qual_to_sanger if $params->{convert};
    return $this;
}

=item convert_qual_to_sanger

Convert the quality string from Illumina to Sanger scaling.

=cut

sub convert_qual_to_sanger {
    my $this=shift;
    return $this->{qual}=join('', map { chr(ord($_)-31) } split(//, $this->{qual}));
}


# PRIVATE NONMEMBER FUNCTION PARSES ILLUMINA READ ID
# ILLUMINA READ ID HAS ':'-SEPARATED FIELDS:
# - UNIQUE INSTRUMENT NAME
# - FLOWCELL LANE
# - TILE NUMBER
# - X-COORDINATE
# - Y-COORDINATE
# ADDITIONALLY, THE FOLLOWING OPTIONAL FIELDS EXIST:
# - #NNNNNNNNNN => BARCODE SEQUENCE
# - /1 OR /2 => FORWARD OR REVERSE PAIRED READ
sub _parse_header {
    my ($this,$hdr,$format)=@_;
    $this->{id}=undef; # complete ID (e.g. "A/1#GATTACA"); always defined
    $this->{base}=undef; # base ID only (e.g. "A"); always defined
    $this->{pair}=undef; # pair ID only (e.g. "1"); only defined if paired
    $this->{barcode}=undef; # barcode sequence (upper-case); only defined it barcoded

    if ($format eq 'illumina') {
        if ($hdr =~ /^@(\S+:\d+:\d+:\d+:\d+)#([aAtTcCgGnN]+)\/([12])/) { # barcoded, paired
            $this->{base}=$1;
            $this->{barcode}=uc($2);
            $this->{pair}=$3;
        } elsif ($hdr =~ /^@(\S+:\d+:\d+:\d+:\d+)\/([12])#([aAtTcCgGnN]+)/) { # barcoded, paired
            $this->{base}=$1;
            $this->{pair}=$2;
            $this->{barcode}=uc($3);
        } elsif ($hdr =~ /^@(\S+:\d+:\d+:\d+:\d+)\/([12])/) { # paired
            $this->{base}=$1;
            $this->{pair}=$2;
        } elsif ($hdr =~ /^@(\S+:\d+:\d+:\d+:\d+)#([aAtTcCgGnN]+)/) { # barcoded, unpaired
            $this->{base}=$1;
            $this->{barcode}=uc($2);
        } elsif ($hdr =~ /^@(\S+:\d+:\d+:\d+:\d+)/) { # unpaired
            $this->{base}=$1;
        } else {
            die("Unable to parse Illumina sequence header: $hdr\n");
        }
        $this->{id}=$this->{base};
        $this->{id} .= "/".$this->{pair} if $this->{pair};
        $this->{id} .= "#".$this->{barcode} if $this->{barcode}
    } elsif ($format eq 'roche') {
        if ($hdr =~ /^@(\S+) length=\d+ xy=\d+_\d+ region=\d run=\S/) {
            $this->{id}=$this->{base}=$1;
        } elsif ($hdr =~ /^@(\S+)/) {
            $this->{id}=$this->{base}=$1;
        } else {
            die("Unable to parse Roche sequence header: $hdr\n");
        }
    } elsif ($hdr =~ /^@(\S+)/) {
        # generic fastq (not an Illumina or Roche read header)
        $this->{id}=$this->{base}=$1;
    } else {
        die("Unable to parse generic sequence header: $hdr\n");
    }
}

# TRIM THE ROCHE MOLECULAR ID FROM THE 5' END OF THE SEQUENCE
sub _trim_roche_mid {
    my ($this,$len)=@_;
    die("Invalid MID length, $len\n") unless $len > 0;
    if (length($this->{seq}) <= $len) {
        $this->{seq}=$this->{qual}=undef;
        return;
    }
    $this->{barcode}=uc(substr($this->{seq},0,$len));
    $this->{seq}=substr($this->{seq},$len);
    $this->{qual}=substr($this->{qual},$len);
    $this->{id}=$this->{base};
    $this->{id} .= "/".$this->{pair} if $this->{pair};
    $this->{id} .= "#".$this->{barcode} if $this->{barcode}
}

# DESTRUCTOR - closes file before going away
sub DESTROY {
    my $this=shift;
    close $this->{fh} if exists($this->{fh}) and defined($this->{fh});
}

=item header ($new_hdr)

Returns the object's header line.  You can also use this method to give the sequence a new header.

=cut

sub header {
    my ($this,$new_hdr)=@_;
    if (defined($new_hdr) and $new_hdr) {
        $new_hdr = $1 if $new_hdr =~ /^>(.+)$/;
        $new_hdr = '@'.$new_hdr unless $new_hdr =~ /^@/;
        $this->_parse_header($new_hdr);
    }
    return '@'.$this->{id};
}

=item id

Returns the object's ID, which is the sequence's unique identifier without comments which may be present in the header.
It cannot be changed directly, but will be updated whenever the header, base, or pair is changed.

=cut

sub id {
    my $this=shift;
    return $this->{id};
}

=item base

If the read is paired, returns it's base ID (same for both members of pair); returns undefined if not paired.
Providing the optional argument will change the base (and the id and header).

=cut

sub base {
    my ($this,$base)=@_;
    if ($base) {
        $this->{id}=$this->{base}=$base;
        $this->{id} .= "/".$this->{pair} if $this->{pair};
        $this->{id} .= "#".$this->{barcode} if $this->{barcode}
    }
    return $this->{base};
}

=item barcode

If the read contains an Illumina molecular barcode ID, it will be returned; otherwise returns undef.
Supplying an optional argument will set the barcode.

=cut

sub barcode {
    my ($this,$barcode)=@_;
    if ($barcode) {
        $this->{barcode}=$barcode;
        $this->{id}=$this->{base};
        $this->{id} .= "/".$this->{pair} if $this->{pair};
        $this->{id} .= "#".$barcode;
    }
    return $this->{barcode};
}

=item pair

Returns the read's ord in the pair; undef otherwise.  It may be changed by supplying the optional extra argument.
To clear the pairing information, pass it an empty string, not NULL.

=cut

sub pair {
    my ($this,$pair)=@_;
    if (defined($pair)) {
        if ($pair) {
            $this->{pair}=$pair;
            $this->{id}=$this->{base};
            $this->{id} .= "/".$pair;
            $this->{id} .= "#".$this->{barcode} if $this->{barcode}
        } else {
            # delete pairing information (e.g. singleton)
            $this->{pair}=undef;
            $this->{id}=$this->{base};
            $this->{id} .= "#".$this->{barcode} if $this->{barcode}
        }
    }
    return $this->{pair};
}

=item seq ($new_seq)

Returns the read's complete sequence, without newlines.  Optional argument changes it.  It's up to the developer to make
sure the quality string is also updated, if necessary!

=cut

sub seq {
    my ($this,$new_seq)=@_;
    $this->{seq}=$new_seq if $new_seq;
    return $this->{seq};
}

=item revcomp

Reverse-complements a sequence and quality scores.

=cut

sub revcomp {
    my $this=shift;
    return unless $this->{seq};
    $this->{seq} =~ tr/ATCGatcg/TAGCtagc/;
    my @seq=reverse split(//, $this->{seq});
    $this->{seq}=join('', @seq);
    my @qual=reverse split(//, $this->{qual});
    $this->{qual}=join('', @qual);
}

=item qual ($new_qual)

Returns the read's quality string, without newlines.  Optional argument changes it.  It's up to the developer to make sure
the seq is also updated, if necessary!

=cut

sub qual {
    my ($this,$new_qual)=@_;
    $this->{qual}=$new_qual if $new_qual;
    return $this->{qual};
}

=item qual_arrayref

Returns an arrayref of sanger-scaled quality scores.

=cut

sub qual_arrayref {
    my $this=shift;
    return [] unless $this->{qual};
    my @qual=map { ord($_)-33 } split(//, $this->{qual});
    return \@qual;
}

=item output

Returns a multiline string of the sequence in Fastq format.  Returns no output if sequence is empty.
Optional convert flag indicates output should have Illumina-scaled quality scores.

=cut

sub output {
    my ($this,$convert)=@_;
    return '' unless $this->{seq}; # will be empty if failed QC
    my $qual = $this->{qual};
    $qual = join('', map { chr(ord($_)+31) } split(//, $this->{qual})) if $convert;
    return '@'.$this->{id}."\n"._format($this->{seq})."+\n"._format($qual);
}

# PRIVATE METHOD TO ADD NEWLINE CHARACTERS 
sub _format {
    my $old=shift;
    return '' unless $old;
    return $old."\n" unless length($old) > CHARACTERS_PER_LINE;
    my $new='';
    while (length($old)> CHARACTERS_PER_LINE) {
        $new .= substr($old,0, CHARACTERS_PER_LINE)."\n";
        $old = substr($old, CHARACTERS_PER_LINE);
    }
    $new .= $old."\n";
    return $new;
}

=item output_fasta

Returns a multiline string of the sequence in Fasta format.  Returns empty string if sequence is empty.

=cut

sub output_fasta {
    my $this=shift;
    return '' unless $this->{seq};
    return '>'.$this->{id}."\n"._format($this->{seq})."\n";
}

=item output_qual

Returns a multiline string of the sequence's quality scores in phred format.  Returns empty string if no scores.

=cut

sub output_qual {
    my $this=shift;
    return '' unless $this->{qual};
    my @qual=map { ord($_)-33 } split(//, $this->{qual});
    my $output=">".$this->{id}."\n";
    my $i= CHARACTERS_PER_LINE - 1;
    while (scalar(@qual) > CHARACTERS_PER_LINE) {
        $output .= join(' ', @qual[0 .. $i])."\n";
        @qual=@qual[CHARACTERS_PER_LINE .. $#qual];
    }
    $output .= join(' ', @qual)."\n";
    return $output
}

###############################################################################
## QC METHODS

=back

=head2 QC Methods

=over 5

=item qc $winsize $meanq $minlen $maxn

To perform quality end trimming, minimum length filtering, and filtering reads with too many Ns.

=cut

sub qc {
    my ($this, $winsize, $meanq, $minlen, $maxn, $low_complexity)=@_;
    $this->trim_terminal_Ns;
    $this->qual_end_trim($winsize, $meanq);
    $this->length_filter($minlen);
    $this->N_filter($maxn);
    $this->low_complexity_filter($low_complexity);
}

=item trim_terminal_Ns

Removes uninformative Ns at the sequence ends.

=cut

sub trim_terminal_Ns {
    my $this=shift;
    if ($this->{seq} =~ /^(N+)(.*)$/) {
        $this->{seq}=$2;
        $this->{qual}=substr($this->{qual},length($1));
    }
    if ($this->{seq} =~ /^(.*)N+$/) {
        $this->{seq}=$1;
        $this->{qual}=substr($this->{qual},0,length($1));
    }
}

=item qual_end_trim

To trim both ends of the read using a sliding window until the mean quality score surpasses the threshold.
Returns length of resultant sequence.  The sequence can be trimmed away to nothing; the output method will not return
anything if the sequence is empty, so you must handle such cases manually, if necessary.  This method will produce 
erroneous output if the fastq object was initialized without being given the convert flag and the Q-scores are actually
Illumina-scaled.

=cut

sub qual_end_trim {
    my ($this, $winsize, $meanq)=@_;
    $winsize=DEFAULT_WINSIZE unless $winsize;
    $meanq=DEFAULT_MEANQ unless $meanq;
    return 0 unless $this->{seq};

    my @qual=map { ord($_)-33 } split(//, $this->{qual});
    if (@qual < $winsize ) {
        # FAIL
        $this->{seq}=$this->{qual}='';
        return 0;
    }

    # SLIDING WINDOW TRIM LEFT
    my $w=$winsize-1;
    my $winsum=0;
    foreach (@qual[0..$w]) { $winsum += $_ }
    while ( @qual>$winsize and ($winsum/$winsize)<$meanq) {
        $this->{seq} = substr($this->{seq},1);
        $this->{qual} = substr($this->{qual},1);
        $q=shift @qual;
        $winsum -= $q;
        $q=$qual[$w];
        $winsum += $q;
    }
    if ( ($winsum/$winsize)<$meanq) {
        # FAIL
        $this->{seq}=$this->{qual}='';
        return 0;
    }

    # SLIDING WINDOW TRIM RIGHT
    $w=@qual-$winsize;
    $winsum=0;
    foreach (@qual[$w..$#qual]) { $winsum += $_ }
    while ( @qual>$winsize and ($winsum/$winsize)<$meanq) {
        $this->{seq} = substr($this->{seq},0,length($this->{seq})-1);
        $this->{qual} = substr($this->{qual},0,length($this->{qual})-1);
        $q=pop @qual;
        $winsum -= $q;
        $w=@qual-$winsize;
        $q=$qual[$w];
        $winsum += $q;
    }
    if ( ($winsum/$winsize)<$meanq) {
        # FAIL
        $this->{seq}=$this->{qual}='';
        return 0;
    }
    return length($this->{seq});
}

=item length_filter

If the sequence is shorter than the minimum length, the sequence and quality string are emptied so they will not be
returned by the output method.  Returns true if sequence was filtered.

=cut

sub length_filter {
    my ($this,$minlen)=@_;
    $minlen=DEFAULT_MINLEN unless defined($minlen);
    if ( length($this->{seq}) >= $minlen) { 
        return 0;
    } else {
        $this->{seq}=$this->{qual}='';
        return 1;
    }
}

=item N_filter

If the sequence contains more than the allowed number of Ns, the sequence and quality string are emptied.
Returns true if sequence was filtered.

=cut

sub N_filter {
    my ($this,$maxn)=@_;
    $maxn=DEFAULT_MAXN unless defined($maxn);

    return 1 unless defined($maxn);
    # COUNT NUMBER OF 'N' BASES IN TRIMMED SEQ
    my $tmpseq=$this->{seq};
    my $n= $tmpseq=~ s/N//gi;
    if ($n <= $maxn) {
        return 0;
    } else {
        # FAIL
        $this->{seq}=$this->{qual}='';
        return 1;
    }
}

=item low_complexity_filter $pct_len

If the sequence is >= $pct_len mono- or di-nucleotide repeats, clears the sequence and quality strings and returns true.

=cut

sub low_complexity_filter {
    my ($this,$pct_len)=@_;
    $pct_len=DEFAULT_LOW_COMPLEXITY unless defined($pct_len);
    my $seq=$this->{seq};    
    my $len = length($seq);
    return 1 unless $len;
    my $filter=0;
    foreach my $nn (qw/AA TT CC GG CA GT CT GA AT CG/) {
        my $n = $seq =~ s/$nn/$nn/gi;
        if ($n*2/$len >= $pct_len) {
            $filter = 1;
            last;
        }
    }
    if ($filter) {
        $this->{seq}=$this->{qual}='';
    }
    return $filter;
}

#=item trim
#
#Given start and end coordinates of adapter/primer sequence, trims sequence and quality strings and returns final length.
#
#=cut
#
#sub trim {
#    my ($this, $start, $end)=@_;
#    return unless defined($start) and defined($end);
#    ($start,$end)=($end,$start) if $end<$start; 
#    my $len=length($this->{seq});
#    #print STDERR "- trim read, ", $this->id, " ($len bp) from $start-$end\n"; # DEBUG
#
#    if ($start <= CLOSE_ENOUGH_TO_END and $len - $end <= CLOSE_ENOUGH_TO_END ) {
#        # fail read
#        $this->{seq}=$this->{qual}='';
#    } elsif ($start <= CLOSE_ENOUGH_TO_END ) {
#        # trim left
#        $this->{seq}=substr($this->{seq},$end);
#    } elsif ($len - $end <= CLOSE_ENOUGH_TO_END ) {
#        # trim right
#        $this->{seq}=substr($this->{seq},0,$start);
#    }
#    return length($this->{seq});
#    #print STDERR "\t+ length after trimming is ", length($this->{seq})," bp\n"; # DEBUG
#}

=back

=head1 BUGS AND LIMITATIONS

Reads must be named in accordance with Illumina naming conventions.

=head1 COPYRIGHT

Copyright (c) 2010 U.S. Department of Energy Joint Genome Institute

All right reserved. This program is free software; you can redistribute it
and/or modify it under the same terms as Perl itself.

=head1 AUTHOR

Edward Kirton <ESKirton@LBL.gov>

=cut

###############################################################################
## FASTQDB OBJECT

=head1 NAME

FastqDb - Simple object for handling Fastq files

=head1 SYNOPSIS

    # simple script to QC reads (outputs to stdout)
    use FastqDb;
    my $file=shift or die("Infile required\n");
    my $fastqFile = new FastqFile($file);
    while (my $fastq=$fastqFile->next_seq) {
        $fastq->qc; # see Fastq package for details
        print $fastq->output; # no output if read failed QC
    }

=head1 DESCRIPTION

This module provides an iterator for Fastq files, plus basic QC functions

=head1 METHODS

=over 5

=cut

package FastqDb;

use IO::File;
use File::Basename;
use Cwd;

=item new $file $convert

Constructor requires a path to an existing Fastq file.

Optional parameters are passed in a hashref:

    convert => if true, convert qual from illumina to sanger encoding
    format => 'roche' or '454' or 'illumina' or 'generic'

If the optional convert parameter is not defined, the first 1000 sequences in the file will be checked in order to infer the
quality encoding used.  Illumina quality scores are re-scaled to the Sanger encoding.

The default format is 'generic', which means barcodes and pairing is not supported.  With Illumina format, the pairing
is encoded by read IDs with "/1" or "/2" and barcodes are encoded in the read ID following the "#" character.  The
formats "roche" and "454" are synonymous; barcodes are located in the 5' end of the read sequence; pairing is not
currently supported.

=cut

sub new {
    my ($class,$file,$convert,$params)=@_;

    # VALIDATE
    die("Fastq file required\n") unless $file;
    $params={ format => 'generic' } unless defined($params);

    # USE FULL PATH
    my $cwd=getcwd.'/';
    my ($base,$dir,$suffix)=fileparse($file, qr/\.[^.]*/);
    $dir=$cwd if $dir eq './';
    $file="$dir$base$suffix";
 
    # INIT
    my $this = { file => $file, convert => undef };
    $this->{fh} = new IO::File("<$file") or die("Unable to open Fastq file, $file: $!\n");
    $this->{params}=$params;
    bless $this,$class;

    my %valid_formats=( 'generic' => undef, 'illumina' => undef, 'roche' => undef, '454' => undef );
    foreach my $param (keys %{$this->{params}}) {
        $param=lc($param);
        $this->{params}->{$param} = lc($this->{params}->{$param});
        if ($param eq 'format') {
            die("Invalid format parameter, $format\n") unless exists($valid_formats{$this->{params}->{format}});
        } elsif ($param eq 'roche_mid_len') {
            my $len=$this->{params}->{roche_mid_len};
            die("Invalid Roche MID length, $len\n") unless $len > 0;
        } else {
            die("Unknown parameter, $param\n");
        }   
    }

    # check qual format if necessary
    if (defined($convert)) {
        if ($convert == 0 or $convert == 1) {
            $this->{params}->{convert}=$convert;
        } else {
            die("Invalid convert paramter, $convert\n");
        }
    } else {
        $this->{params}->{convert}= $this->is_sanger ? 0:1;
    }

    return $this;
}

=item is_sanger

This method will return true only if the file has sanger-encoded quality scores, but it will not set the object's
'convert' flag or alter the file handle's current position.

=cut

sub is_sanger {
    my $this=shift;
    # Check first 1000 records to determine quality score encoding method (Illumina 1.3+ or Sanger; Solexa not supported)
    my $counter=0;
    my $appears_sanger=0;
    my $appears_illumina=0;
    my $fh=$this->{fh};
    $this->{fh} = new IO::File('<'.$this->{file});
    while (my $seq=$this->next_seq) {
        my @qual=split(//, $seq->qual);
        foreach my $q (@qual) {
            my $ascii=ord($q);
            if ( $ascii >= 35 and $ascii <= 65 ) { ++$appears_sanger }
            elsif ( $ascii >= 74 and $ascii <= 104 ) { ++$appears_illumina }
        }
        last if ++$counter > 1000;
    }
    $this->close_file;
    $this->{fh}=$fh;
    return $appears_sanger > $appears_illumina ? 1:0;
}

=item file ($new_file)

Returns the Fastq db complete path.  You may also use this method to open a new file for reading.

=cut

sub file {
    my ($this,$new_file)=@_;
    if ($new_file) {
        $this->close_file;
        $this->{fh}=new IO::File('<'.$new_file) or die("Unable to open Fastq file, $new_file: $!\n");
    }
    return $this->{file};
}

=item convert ($flag)

Returns the convert flag's value; 0=do not convert; 1=convert from Illumina to Sanger-scaling.  
Optionally, you may use this method to set the flag although the check_qual_format method is preferred;
setting the convert flag for a sanger-encoded file will result in nonsense output!

=cut

sub convert {
    my ($this,$convert)=@_;
    if (defined($convert)) { $this->{convert} = $convert ? 1:0 }
    return $this->{convert};
}

=item next_seq $convert_to_illumina

An iterator for the fastq database; returns a single sequence as multiline text.

=cut

sub next_seq {
    my $this=shift;
    return undef unless exists($this->{fh}) and defined($this->{fh});

    # INIT VARS
    my $fh = $this->{fh};
    my $hdr=undef; # complete header line without \n
    my $seq=''; # without \n
    my $qual=''; # without \n

    # GET ID
    while (<$fh>) {
        chomp;
        if (/^#/ or ! $_) {
            next;
        } elsif (/^@\S+/) {
            $hdr=$_;
        } else {
            die("Expected Fastq header line, got: \"$_\"\n");
        }
        last;
    }
    unless ($hdr) {
        $this->close_file;
        return undef;
    }

    # GET SEQ
    while (<$fh>) {
        chomp;
        if (/^\+/) { last } else { $seq .= $_ }
    }
    die("Missing seq for $hdr\n") unless $seq;

    # GET QUAL
    my $ls=length($seq);
    my $lq=0;
    while (<$fh>) {
        chomp;
        $qual .= $_;
        $lq += length($_);
        if ($lq == $ls) {   
            last;
        } elsif ($lq > $ls) {
            die("Invalid record for $hdr; qual string is longer than sequence string:\n\t$seq\n\t$qual\n");
        }
    }
    die("Incomplete record for $hdr\n") unless $qual;

    # RETURN FASTQ OBJECT
    return new Fastq( $hdr, $seq, $qual, $this->{params} );
}

=item close_file

This will close the file and delete the file handle.  Normally this is called by the iterator upon reaching the end
of the file, so this method only needs to be called if you wish to close the file prematurely.

=cut

sub close_file {
    my $this=shift;
    return unless exists($this->{fh}) and defined($this->{fh});
    close $this->{fh};
    delete $this->{fh};

}

1;

=back

=head1 BUGS AND LIMITATIONS

Expects read IDs to follow Illumina's naming convention.

=head1 COPYRIGHT

Copyright (c) 2010 U.S. Department of Energy Joint Genome Institute

All right reserved. This program is free software; you can redistribute it
and/or modify it under the same terms as Perl itself.

=head1 AUTHOR

Edward Kirton <ESKirton@LBL.gov>

=cut


