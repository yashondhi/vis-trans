#!/usr/bin/perl

my $input = $ARGV[0];
my $format = $ARGV[1];
my $missing = $ARGV[2];
my $output = "output";
my $fparam;

if($missing eq 'yes'){
	$fparam = "-F";
}else{
	$fparam = "";
}
my $run = qx/prank -d=$input -o=$output -f=$format $fparam/;
print $run;
