#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';

## Input argument : directory path to install
my $in_dir = $ARGV[0];
my $options = GetOptions(
				"dir=s" => \$in_dir
);

#---------------------------------------------------------------------
###### Set up directory
if(!defined($in_dir)){
		$in_dir = abs_path("$Bin/../bin");
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
