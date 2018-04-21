
#!/usr/bin/perl

my $new_sequences = $ARGV[0];
my $existing_alignment = $ARGV[1];
my $output = "seqs_aligned.fasta";

system "mafft --add $new_sequences --reorder $existing_alignment > $output 2>log.txt ";

