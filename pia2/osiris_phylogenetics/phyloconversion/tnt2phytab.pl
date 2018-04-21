#!/usr/bin/perl

use strict;

my $file = $ARGV[0];

# read file with genes
open FILE, $file or die "ERROR: Cannot open file $file\n";
my $firstline=0;
my $datatype;
my $taxa;
while (<FILE>) {
       my $currentinput = "$_";
       if($firstline==0){
               if($currentinput =~ m/nstates/){
                       my @splitlines=split(' ',$currentinput);
                       $splitlines[2] =~ s/\;//;
                       if($splitlines[2] == 2){
                               $datatype = "binary";
                       }elsif($splitlines[2] > 2){
                               $datatype = "multi";
                       }
               }else{
                       die "ERROR: file does not begin with nstates line. Must be TNT file exported from MorphoBank.org";
               }
       }
       if($firstline==1){
               if($currentinput =~ m/xread/){
               }else{
                       die "ERROR: file does not contain xread line. Must be TNT file exported from MorphoBank.org";
               }
       }
       if($firstline==2){
               if($currentinput =~ m/Morpho/){
               }else{
                       die "ERROR: file does not contain Morphobank Comment line. Must be TNT file exported from MorphoBank.org";
               }
       }
       if($firstline==3){
               if($currentinput =~ m/\d/){
                       my @splitlines=split(' ',$currentinput);
                       $taxa = $splitlines[1]."\n";
               }else{
                       die "ERROR: file does not contain number of taxa. Must be TNT file exported from MorphoBank.org";
               }
       }
       if($firstline==4){
               if($currentinput =~ m/\d/){
                       die "ERROR: file does not contain empty line after taxa numbers . Must be TNT file exported from MorphoBank.org";
               }else{
               }
       }
       if(($firstline>4)&&($firstline<(3+2+$taxa))){
               my @splitlines=split(' ',$currentinput);
               print $splitlines[0]."\t".$datatype."\t".$splitlines[0]."_".$datatype."\t".$splitlines[1]."\n";
       }
       $firstline++;
}
