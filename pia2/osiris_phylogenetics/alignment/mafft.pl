#!/usr/bin/perl

my $strategy = $ARGV[0];
my $input = $ARGV[1];
my $output = "seqs_aligned.fasta";

if($strategy eq "Auto") {
	my $run = qx/mafft --auto $input > $output 2>log.txt/;
}
elsif($strategy eq "FFT-NS-1") {
	my $run = qx/mafft --retree 1 $input > $output 2>log.txt/;
}
elsif($strategy eq "FFT-NS-2") {
	my $run = qx/mafft --retree 2 $input > $output 2>log.txt/;
}
elsif($strategy eq "FFT-NS-i") {
	my $run = qx/mafft-fftnsi $input > $output 2>log.txt/;
}
elsif($strategy eq "E-INS-i") {
	my $run = qx/mafft-einsi $input > $output 2>log.txt/;
}
elsif($strategy eq "L-INS-i") {
	my $run = qx/mafft-linsi $input > $output 2>log.txt/;
}
elsif($strategy eq "G-INS-i") {
	my $run = qx/mafft-ginsi $input > $output 2>log.txt/;
}
elsif($strategy eq "Q-INS-i") {
	my $run = qx/mafft-qinsi $input > $output 2>log.txt/;
}

print $run;
