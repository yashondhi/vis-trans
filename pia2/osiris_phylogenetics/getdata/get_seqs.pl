#!/jgi/tools/bin/perl -w

#
# This script creates a Fasta/Qual/Fastq file of selected sequences, with optional filters.
#
# 02/24/10 : created by Ed Kirton
# 12/07/10 : fixed Fastq bug
#

use strict;
use warnings;
use Getopt::Long;
use IO::File;
#use PerlIO::gzip;
use FindBin;
use lib $FindBin::Bin;
use FastaDb;
use FastqDb;

my $usage = <<'ENDHERE';
NAME:
    get_seqs.pl
PURPOSE:
    To extract a subset of sequences by ID.
INPUT:
    --db <*.fasta|fastq> : file containing sequences in Fasta or Fastq format
    --table <*.tsv> : file containing sequence IDs (optional; default=stdin)
    --col <int> : column of table containing sequence IDs (optional; default=1=first column)
OUTPUT:
    --selected <*.fasta|fastq> : file containing named sequences
    --unselected <*.fasta|fastq> : file containing unselected sequences
OPTIONS:
    --cosorted : uses faster algorithm if IDs appear in both files in the same order
    --paired : filter complete read-pair when one read is selected (requires Illumina-style read IDs; i.e. */1, */2)
    --ignore case : ignore case differences between IDs
    --gzip : compress all outfiles
OPTIONAL FILTERS:
    Optional filters, each in the form of column:condition:value.
    Where column is the column in the table (containing IDs)
    Condition is one of the following:
        String operators:
            s_eq
            s_ne
            s_contains
            s_notcontains
            s_startswith
            s_notstartswith
            s_endswith
            s_notendswith
        Numerical operators:
            n_eq
            n_ne
            n_gt
            n_lt
    Where value is a string or number as appropriate.
AUTHOR/SUPPORT:
    Edward Kirton (ESKirton@LBL.gov)
ENDHERE

#
# VALIDATE INPUT
#
my ($help, $dbfile, $tablefile, $id_col, $ignorecase, $cosorted, $selected, $unselected, $gzip, $paired);
GetOptions(
    'd|db=s'           => \$dbfile,
    't|table=s'        => \$tablefile,
    'c|col=i'          => \$id_col,
    'ignorecase'     => \$ignorecase,
    'cosorted'       => \$cosorted,
    's|selected=s'   => \$selected,
    'u|unselected=s' => \$unselected,
    'g|gzip' => \$gzip,
    'p|paired' => \$paired,
    'h|help'         => \$help
);
if ($help) { print $usage; exit; }
die("DB required\n") unless $dbfile;
die("DB file not found: $dbfile\n") unless -f $dbfile;
die("Table required\n") unless $tablefile;
die("Table file not found: $tablefile\n") unless -f $tablefile;
$selected   = '' if !defined($selected) or $selected   eq 'None';
$unselected = '' if !defined($unselected) or $unselected eq 'None';
$id_col=1 unless $id_col;
die("Invalid id column, $id_col\n") unless $id_col > 0;

my $filters = [];
while (my $filter = shift @ARGV) {
    next unless $filter;
    my @a_filter = split(/:/, $filter);
    die("Invalid number of filter options: @a_filter") unless @a_filter == 3;
    push @$filters, \@a_filter;
}

#
# MAIN
#
my ($n_selected,$n_unselected);
if ($cosorted) {
    # SEARCH IS FAST AND EASY IF INPUTS SIMILARLY SORTED!
    ($n_selected,$n_unselected) = search_cosorted($dbfile, $tablefile, $id_col, $ignorecase, $selected, $unselected, $paired, $gzip, $filters);
} else {
    # INPUT NOT CO-SORTED SO KEEP ALL IDS IN RAM
    ($n_selected,$n_unselected) = search($dbfile, $tablefile, $id_col, $ignorecase, $selected, $unselected, $paired, $gzip, $filters);
}
print "Selected = $n_selected; Unselected = $n_unselected\n"; 
exit;

#
# RETURNS TRUE ONLY IF RECORD MATCHES (OPTIONAL) SEARCH CRITERIA
#
sub match
{
    my ($filters, $row) = @_;
    foreach my $filterA (@$filters) {
        my ($condition, $col, $value) = @$filterA;
        my $x = $row->[ $col - 1 ];
        if ($condition eq 's_eq') { return 0 unless $x eq $value }
        elsif ($condition eq 's_ne') { return 0 unless $x ne $value }
        elsif ($condition eq 's_contains') { return 0 unless $x =~ /$value/ }
        elsif ($condition eq 's_notcontains')   { return 0 unless $x !~ /$value/ }
        elsif ($condition eq 's_startswith')    { return 0 unless $x =~ /^$value/ }
        elsif ($condition eq 's_notstartswith') { return 0 unless $x !~ /^$value/ }
        elsif ($condition eq 's_endswith')      { return 0 unless $x =~ /$value$/ }
        elsif ($condition eq 's_notendswith')   { return 0 unless $x !~ /$value$/ }
        elsif ($condition eq 'n_eq')            { return 0 unless $x == $value }
        elsif ($condition eq 'n_ne')            { return 0 unless $x != $value }
        elsif ($condition eq 'n_gt')            { return 0 unless $x > $value }
        elsif ($condition eq 'n_lt')            { return 0 unless $x < $value }
    }
    return 1;
}

#
# SIMULTANEOUSLY PARSE TWO STREAMS
#
sub search_cosorted
{
    my ($dbfile, $tablefile, $id_col, $ignorecase, $selected, $unselected, $paired, $gzip, $filters) = @_;
    my $sfh = new IO::File;
    my $ufh = new IO::File;
    my $table = new IO::File;
    my $n_selected = 0;
    my $n_unselected = 0;

    # OPEN FILES
    if ($tablefile) {
        open($table, "<$tablefile") or die("Unable to open file, $tablefile: $!\n");
    } else {
        $table=*STDIN;
    }
    if ($selected) {
        if ($gzip) {
#            open($sfh, '>:gzip', $selected) or die("Unable to open file, $selected: $!\n");
        } else {
            open($sfh, ">$selected") or die("Unable to open file, $selected: $!\n");
        }
    } else {
        open($sfh, ">/dev/null");
    }
    if ($unselected) {
        if ($gzip) {
#            open($ufh, '>:gzip', $unselected) or die("Unable to open file, $unselected: $!\n");
        } else {
            open($ufh, ">$unselected") or die("Unable to open file, $unselected: $!\n");
        }
    } else {
        open($ufh, ">/dev/null");
    }

    # GET FIRST MATCHING TARGET ID
    my $prev_target_id = '';
    my $target_id = '';
    get_next_matching_target_id($table,$id_col,$ignorecase,$filters,\$target_id,\$prev_target_id,$paired);
    unless ($target_id) {
        # no records match search criteria
        close $table;
        close $sfh if $selected;
        if ($unselected) {
            open(DB, "<$dbfile") or die("Unable to open file, $dbfile: $!\n");
            while (<DB>) {
                print $ufh $_;
                ++$n_unselected;
            }
            close DB;
        }
        close $ufh;
        return 0;
    }

    # DETERMINE FILETYPE
    open(DB, "<$dbfile") or die("Unable to open file, $dbfile: $!\n");
    my $format;
    while (<DB>) {
        chomp;
        if (/^#/ or ! $_) { next }
        elsif (/^>/) { $format='fasta' }
        elsif (/^@/) { $format='fastq' }
        else { die "Invalid DB file format" }
        last;
    }
    close DB;

    # PARSE
    my $db = $format eq 'fasta' ? FastaDb->new($dbfile) : FastqDb->new($dbfile);
    while (my $rec=$db->next_seq ) {
        unless ($target_id) {
            last unless $unselected; # done if no more seqs to get
            # otherwise dump rest of seqs in unselected file
            print $ufh $rec->output;
            ++$n_unselected;
            while ($rec=$db->next_seq ) {
                print $ufh $rec->output;
                ++$n_unselected;
            }
            last;
        }
        my $id=$ignorecase ? uc($rec->id):$rec->id;
        if ($id eq $prev_target_id or $id eq $target_id) {
            # selected seq
            print $sfh $rec->output;
            ++$n_selected;
            get_next_matching_target_id($table,$id_col,$ignorecase,$filters,\$target_id,\$prev_target_id,$paired);
        } else {
            # unselected seq
            print $ufh $rec->output;
            ++$n_unselected;
        }
    }
    close $table;
    close $sfh;
    close $ufh;

    # If some target seqs not found, it's likely the files were not cosorted, so try unsorted search function.
    if ($target_id) {
        print "Files don't appear to be cosorted, trying unsorted search\n";
        return search($dbfile, $tablefile, $id_col, $ignorecase, $selected, $unselected, $filters);
    }
    return ($n_selected,$n_unselected);

    sub get_next_matching_target_id {
        my ($table,$id_col,$ignorecase,$filters,$target_idR,$prev_target_idR,$paired)=@_;
        $$prev_target_idR = $$target_idR;
        $$target_idR = '';
        while (<$table>) {
            chomp;
            my @row = split(/\t/);
            die("Bad input table") unless @row >= $id_col;
            next unless match($filters, \@row);
            my $new_target_id = $ignorecase ? uc($row[ $id_col - 1 ]) : $row[ $id_col - 1 ];
            $new_target_id=$1 if $new_target_id =~ /^(\S+)/; # use first word only
            $new_target_id=$1 if $paired and $new_target_id =~ /^(\S+)\/[12]$/;
            next if $new_target_id eq $$prev_target_idR;
            $$target_idR=$new_target_id;
            last;    # return to parsing db file
        }
    }
}

#
# LOAD IDS IN RAM THEN PARSE DB.
#
sub search
{
    my ($dbfile, $tablefile, $id_col, $ignorecase, $selected, $unselected, $paired, $gzip, $filters) = @_;
    my $sfh = new IO::File;    # selected seqs
    my $ufh = new IO::File;    # unselected seqs
    my $table=new IO::File;
    my $n_selected=0;
    my $n_unselected=0;
    my %ids = ();
    open(DB,    "<$dbfile")    or die("Unable to open file, $dbfile: $!\n");
    if ($tablefile) {
        open($table, "<$tablefile") or die("Unable to open file, $tablefile: $!\n");
    } else {
        $table=*STDIN;
    }
    if ($selected) {
        if ($gzip) {
#            open($sfh, '>:gzip', $selected) or die("Unable to open file, $selected: $!\n");
        } else {
            open($sfh, ">$selected") or die("Unable to open file, $selected: $!\n");
        }
    } else {
        open($sfh, ">/dev/null");
    }
    if ($unselected) {
        if ($gzip) {
#            open($ufh, '>:gzip', $unselected) or die("Unable to open file, $unselected: $!\n");
        } else {
            open($ufh, ">$unselected") or die("Unable to open file, $unselected: $!\n");
        }
    } else {
        open($ufh, ">/dev/null");
    }

    # LOAD IDS OF MATCHING ROWS
    my $num_targets=0;
    while (<$table>) {
        next if /^#/;
        chomp;
        my @row = split(/\t/);
        my $id = $ignorecase ? uc($row[ $id_col - 1 ]) : $row[ $id_col - 1 ];
        $id=$1 if $id =~ /^(\S+)/;
        $id=$1 if $paired and $id =~ /^(\S+)\/[12]$/;
        if (match($filters, \@row)) {
            # remember this ID
            $ids{$id} = 0; # number of reads with this ID found (counter for paired option)
            ++$num_targets;
        }
    }
    unless ($num_targets) {
        # no records match search criteria
        close $table;
        close $sfh if $selected;
        if ($unselected) {
            open(DB, "<$dbfile") or die("Unable to open file, $dbfile: $!\n");
            while (<DB>) {
                print $ufh $_;
                ++$n_unselected;
            }
            close DB;
        }
        close $ufh;
        return 0;
    }


    # DETERMINE FILETYPE
    open(DB, "<$dbfile") or die("Unable to open file, $dbfile: $!\n");
    my $format;
    while (<DB>) {
        chomp;
        if (/^#/ or /^$/) { next }
        elsif (/^>/) { $format='fasta' }
        elsif (/^@/) { $format='fastq' }
        else { die "Invalid DB file format" }
        last;
    }
    close DB;

    # GET SEQS
    my $db = $format eq 'fasta' ? FastaDb->new($dbfile) : FastqDb->new($dbfile);
    while (my $rec=$db->next_seq ) {
        my $id = $ignorecase ? uc($rec->id) : $rec->id;
        $id = $1 if $paired and $id =~ /^(\S+)\/[12]$/;
        if (exists($ids{$id})) {
            # selected
            print $sfh $rec->output;
            ++$n_selected;
            if (!$paired) {
                delete $ids{$id};
            } else {
                $ids{$id} += 1;
                delete $ids{$id} if $ids{$id} == 2;
            }
        } else {
            # unselected
            print $ufh $rec->output;
            ++$n_unselected;
        }
    }
    close $table;
    close $sfh;
    close $ufh;

    # MAKE SURE ALL TARGETS WERE FOUND
    foreach my $id (keys %ids) {
        if ($ids{$id}) {
            delete($ids{$id}); # SOMETIMES INFILES CONTAIN ONLY ONE READ OF PAIR
	}elsif($id eq 'EMPTY') {		#ADDEDBY THO TO ALLOW EMPTY blast results 
						#in workflow using checkempty.pl
	    exit;		
        } else {
            warn("Seq not found: $id\n");
        }
    }
    return ($n_selected,$n_unselected);
}

__END__
