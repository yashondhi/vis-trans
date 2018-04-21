#!/usr/bin/perl

# UCSB Hamster Installation Script

use strict;
use warnings;
use Cwd;

my $galaxy_hamster = getcwd();

my @tokens = split('/', $galaxy_hamster);

# Checking if all dependencies in environment variable or not.
# If conflicts are found, install will abort. Otherwise, the 
# hamster dependencies will be added to the path

chdir("\/$tokens[1]\/$tokens[2]\/");

open(BASHRCIN, "<".".bashrc");
	while(my $currLine = <BASHRCIN>) {
		if($currLine eq  "export PATH=\"$galaxy_hamster:\$PATH\";\n") {
			die "Error: Please make sure there are no conflicting files in your environment variables of your .bashrc file.\n\nIf you are unsure what to do, refer to the detailed installation procedures in the README file.\nAborted installation.";
		}
	}
close(BASHRCIN);

open(BASHRC, ">>".".bashrc");
	print BASHRC "export PATH=\"$galaxy_hamster:\$PATH\";\n";
close(BASHRC);

# edit wisecfg to include $galaxy_hamster/lib/wisecfg
chdir("$galaxy_hamster\/lib");

open(WISECFG, "<"."run_genewise.pm");
	my $new_run_genewise_pm = "";
	my $i = 0;
	while(my $currLine = <WISECFG>) {
		$i++;
		if($i == 3) {
			$currLine = "$ENV{'WISECONFIGDIR'} = "."\"$galaxy_hamster".'/lib/wisecfg"'.";\n";
		}
		$new_run_genewise_pm = $new_run_genewise_pm.$currLine;
	}
close(WISECFG);

open(WISECFGOUT, ">"."run_genewise.pm");
	print WISECFGOUT $new_run_genewise_pm;
close(WISECFGOUT);

# installation was a success
print "Installation was a SUCCESS!\n\nYou may re-run this script in the future to reconfigure UCSB Hamster.\n\n"
