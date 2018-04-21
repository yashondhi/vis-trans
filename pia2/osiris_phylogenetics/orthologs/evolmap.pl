#!/usr/bin/perl

my $evolmapPath = "\/home\/galaxy\/galaxy_dist\/tool-data\/shared\/jars\/evolmap";

my $file = "Options.txt";	     

my $tree_Input = $ARGV[0];
my $tag_Input = ".dat";
my $protein_Input = $ARGV[1];
my $database_name_Input = "dataout";
my $read_database_Input = $ARGV[3];
my $Blastall_Input = $ARGV[4];
my $read_blast_scores_Input = $ARGV[5];
my $alignments_Input = $ARGV[6];
my $bit_scores_Input = $ARGV[7];
my $read_scores_Input = $ARGV[8];
my $read_ancestors_Input = $ARGV[9];
my $view_ancestors_Input = "false";
my $sfa_Input = $ARGV[10];
my $ortholog_threshold_Input = $ARGV[11];
my $diverged_threshold_Input = $ARGV[12];
my $diverged_std_Input = $ARGV[13];
my $avg_of_paralogs_Input = $ARGV[14];
#my $numDiffGenes = $ARGV[15];

my $temp = $tree_Input;
$temp =~ tr/(),/ /;
my @genes = split(' ', $temp);
my $size = @genes;

my @treeFiles;
my $argIndex = 15;
my $count;
for($count = 0; $count < $size; $count++) {
	$treeFiles[$count] = $ARGV[$argIndex];
	$argIndex++;
}

my $tree_copy = "";
my $index = 0;
my $flag = 1;
for($count = 0; $count < length($tree_Input); $count++) {
	if(substr($tree_Input, $count, 1) eq '(' || substr($tree_Input, $count, 1) eq ')' || substr($tree_Input, $count, 1) eq ',') {
		$tree_copy = $tree_copy.substr($tree_Input, $count, 1);
		$flag = 1;
	}
	else {
		if($flag) {
			$tree_copy = $tree_copy.$treeFiles[$index];
			$index++;
		}
		$flag = 0;
	}
}

$tree_copy =~ s/\Q.dat\E//g;

open(CONFIG, '>'.$file);

print CONFIG "processors = 10\n";
#print CONFIG "tree = ".$tree_Input."\n";
print CONFIG "tree = ".$tree_copy."\n";
print CONFIG "tag = ".$tag_Input."\n";
print CONFIG "protein = ".$protein_Input."\n";
print CONFIG "database_name = ".$database_name_Input."\n";
print CONFIG "read_database = ".$read_database_Input."\n";
print CONFIG "Blastall = ".$Blastall_Input."\n";
print CONFIG "read_blast_scores = ".$read_blast_scores_Input."\n";
print CONFIG "alignments = ".$alignments_Input."\n";
print CONFIG "bit_scores = ".$bit_scores_Input."\n";
print CONFIG "read_scores = ".$read_scores_Input."\n";
print CONFIG "read_ancestors = ".$read_ancestors_Input."\n";
print CONFIG "view_ancestors = ".$view_ancestors_Input."\n";
print CONFIG "sfa = ".$sfa_Input."\n";
print CONFIG "ortholog_threshold = ".$ortholog_threshold_Input."\n";
print CONFIG "diverged_threshold = ".$diverged_threshold_Input."\n";
print CONFIG "diverged_std = ".$diverged_std_Input."\n";
print CONFIG "avg_of_paralogs = ".$avg_of_paralogs_Input."\n";

close(CONFIG);

if($read_ancestors_Input eq "false" || $read_scores_Input eq "false") {
	my $run = qx/java -jar -Xms8000m -Xmx8000m $evolmapPath\/EvolMAP.jar $file/;
	print $run;
}

