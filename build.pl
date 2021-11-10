#!/usr/bin/env perl

use strict;
use warnings;
use FindBin '$Bin';
use File::Basename;
use Cwd 'abs_path';

my $mode = shift;

if(! $mode){
	my $src = basename($0);
	print STDERR "\$./$src install\n";
	print STDERR "\$./$src uninstall\n";

	exit;
}

my $src_path = abs_path("$Bin/sources");
my $thirdparty_path = abs_path("$Bin/third_party");
my $log_dir = "$Bin/third_party/logs";
`mkdir -p $log_dir`;

if($mode eq "install"){
	`chmod +x $Bin/example/*.sh`;
	`mkdir -p $Bin/bin`;
	# Unzip bowtie2
	print STDERR ">> Unzip bowtie2...";
	`wget https://github.com/BenLangmead/bowtie2/releases/download/v2.4.2/bowtie2-2.4.2-linux-x86_64.zip -P $src_path/`;
	`unzip $src_path/bowtie2-2.4.2-linux-x86_64.zip -d $thirdparty_path`;
	`cp $thirdparty_path/bowtie2-2.4.2-linux-x86_64/bowtie2* $Bin/bin/`;
	`rm -f $src_path/bowtie2-2.4.2-linux-x86_64.zip`;
	print STDERR "Done\n";

	# Preparing bedtools
	print STDERR ">> Preparing bedtools...";
	`tar xf $src_path/bedtools_2.17.0.orig.tar.gz -C $thirdparty_path`;
	`mv $thirdparty_path/bedtools-2.17.0 $thirdparty_path/bedtools2`;
	`make -C $thirdparty_path/bedtools2 2> $log_dir/bedtools2.log`;
	if(-f "$thirdparty_path/bedtools2/bin/bedtools"){
		`cp $thirdparty_path/bedtools2/bin/bedtools $Bin/bin/`;
		print STDERR "Done\n";
	}else{
		print STDERR "Error\n";
		exit(1);
	}
		# Preparing samtools
	print STDERR ">> Preparing samtools...";
	`tar xf $src_path/samtools-1.9.tar.gz -C $thirdparty_path`;
	chdir("$thirdparty_path/samtools-1.9");
	`configure 2> ../logs/samtools.config.log`;
	`make 2> ../logs/samtools.log`;
	chdir($Bin);
	if(-f "$thirdparty_path/samtools-1.9/samtools"){
			`cp $thirdparty_path/samtools-1.9/samtools $Bin/bin/`;
		print STDERR "Done\n";
	} else {
		print STDERR "Error\n";
		exit(1);
	}

}else{
	chdir($thirdparty_path);
	`rm -rf *`;
	chdir($Bin);
}
