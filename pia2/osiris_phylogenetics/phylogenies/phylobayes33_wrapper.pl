#!/usr/bin/env perl

my $phylobayes_path = '/home/galaxy/pkgs/phylobayes3.3b/exe_lin64/pb';
my $readPB_path = '/home/galaxy/pkgs/phylobayes3.3b/exe_lin64/readpb';

my $fileName = $ARGV[0];
my $nchainInput = $ARGV[1];
my $cycle_bp_trace_comp = $ARGV[2];
my $discrepancies_threshold = $ARGV[3];
my $effective_size_floor = $ARGV[4];
my $jobName = "dataset";

my $burnin = $ARGV[5];
my $sampleInterval = $ARGV[6];

my $run1 = qx/$phylobayes_path -d $fileName -nchain $nchainInput $cycle_bp_trace_comp $discrepancies_threshold $effective_size_floor $jobName 2>errorlog/;
print $run1;

my $list = qx/ls -l/;
print $list;

my $run2 = qx/$readPB_path -x $burnin $sampleInterval $jobName 2>errorlog/;
print $run2;
