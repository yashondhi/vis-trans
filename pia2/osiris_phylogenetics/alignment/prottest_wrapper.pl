#!/usr/bin/perl

use warnings;
use strict;
use Cwd;

my $dir=getcwd();

#protest directory placed in main user path.  Also, changed runProttest
#script to include full path of jar file
my $prottestPath='/home/galaxy/pkgs/ProtTest2.4';

my $input=$ARGV[1];
my $output=$ARGV[3];

system "$prottestPath/runProtTest -i $input -o $output" ;
