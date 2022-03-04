#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';

## Input argument : directory path to install
my $in_dir;
my $help;
my $bin_dir = abs_path("$Bin/../bin");
my $options = GetOptions(
				"dir=s" => \$in_dir,
				"h|help" => \$help
);

if(defined($help)){ PRINT_HELP(); }
#---------------------------------------------------------------------
###### Set up directory
if(!defined($in_dir)){
		print STDERR "-+- Do you want to install TAMA in to $bin_dir ? [Y/N]\n";
		my $rep = <STDIN>;
		if($rep eq "y" | $rep eq "Y"){
				$in_dir = abs_path("$Bin/../bin");
		}
		else{
				print STDERR "-+- Please set up the place to install TAMA.\n";
				PRINT_HELP();
		}
}
else{
		if(!-e $in_dir){ `mkdir -p $in_dir`; }
				$in_dir = abs_path("$in_dir");
}
print STDERR "TAMA will be downloaded at [$in_dir/TAMA]\n";
#---------------------------------------------------------------------
###### Download and intall TAMA
`git clone https://github.com/jkimlab/TAMA.git $in_dir/TAMA`;
chdir("$in_dir/TAMA");
`./setup.pl --install`;

#---------------------------------------------------------------------
###### Download and set up ready-made species-level databases
`./setup.pl --db --species`;
chdir("$Bin/../bin");
`ln -s $in_dir/TAMA ./`;

####### Print help page
sub PRINT_HELP{
		print STDERR "Usage: $0 -dir <directory to install>\n";
		print STDERR "\tIf -dir is not defined, TAMA will be installed at $bin_dir\n";
		exit();
}
