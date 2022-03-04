#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';
use Parallel::ForkManager;

## Mandatory
my ($read_s, $read1, $read2);
my $in_ref_list;

## Optional
my $out_dir = "./ref_dir";
my $in_cpu = 1;
my $help;

## parameters
my $options = GetOptions(
		"s=s" => \$read_s, # mandatory (incompatible with -1 and -2)
		"1=s" => \$read1, # mandatory (incompatible with -s)
		"2=s" => \$read2, # mandatory (incompatible with -s)
		"o=s" => \$out_dir, # default: ./PR.out
		"p|core=i" => \$in_cpu, # default: 1
		"h|help" => \$help # print help
);

## Parameter check
if(!defined($options)){ HELP(); }
if(defined($help)){ HELP(); }
if(defined($read_s) && defined($read1)){
		HELP();
}
if(defined($read1) && !defined($read2)){
		HELP();
}
if(defined($read2) && !defined($read1)){
		HELP();
}

`mkdir -p $out_dir`;

#---------------------------------------------------------------------
###### Creating a TAMA input parameter file
open(FPR, ">$out_dir/params");
PRINT_PARAM1();
if(defined($read_s)){
		$read_s = abs_path($read_s);
		print FPR "\$SINGLE=$read_s\n\n";
}
else{
		$read1 = abs_path($read1);
		$read2 = abs_path($read2);
		print FPR "\$PAIRED1=$read1\n";
		print FPR "\$PAIRED2=$read2\n\n";
}
PRINT_PARAM2();
close(FPR);

#---------------------------------------------------------------------
###### Running TAMA
`$Bin/../bin/TAMA/TAMA.pl -p $in_cpu -o $out_dir/tamaout --param $out_dir/params`;
my %hs_taxid = ();
#---------------------------------------------------------------------
###### Storing taxonomy IDs of predicted species
open(FAB, "$out_dir/tamaout/TEST/sample1/abundance_profile.0.34.out");
while(<FAB>){
		chomp;
		if($_ =~ /^Scientific/){ next; }
		if($_ =~ /^NA/){ next; }
		my @t = split(/\t/,$_);
		$hs_taxid{$t[1]} = 1;
}
close(FAB);
#---------------------------------------------------------------------
###### Downloading reference genomes with "Chromosome" and "Complete Genome" level from NCBI
`mkdir -p $out_dir/ref_fa`;
open(FTXT, "gunzip -c $Bin/../sources/assembly_summary.txt.gz |");
while(<FTXT>){
		chomp;
		if($_ =~ /^#/){ next; }
		my @t = split(/\t/,$_);
		if($t[19] eq "na"){ next; }
		if(exists($hs_taxid{$t[6]})){
				if($t[11] eq "Chromosome" ||  $t[11] eq "Complete Genome"){
						my @base = split(/\//,$t[19]);
						print STDERR "wget -c -P $out_dir/ref_fa $t[19]/$base[-1]_genomic.fna.gz\n";
						`wget -c -P $out_dir/ref_fa $t[19]/$base[-1]_genomic.fna.gz`;
				}
		}
}
close(FTXT);
#---------------------------------------------------------------------
###### Finishing
`ls $out_dir/ref_fa/* > $out_dir/ref_list.txt`;
#---------------------------------------------------------------------

###### Printing the first part of the parameter file of TAMA
sub PRINT_PARAM1{
		print FPR "[Project]\n";
		print FPR "\$PROJECTNAME=TEST\n\n";
		print FPR "[Basic_options]\n";
		print FPR "\$TOOL=CLARK,centrifuge,kraken\n";
		print FPR "\$RANK=species\n";
		print FPR "\$WEIGHT-CLARK=1\n";
		print FPR "\$WEIGHT-centrifuge=1\n";
		print FPR "\$WEIGHT-kraken=1\n";
		print FPR "\$META-THRESHOLD=\n\n";
		print FPR "[Database]\n";
		print FPR "\$DBNAME=tama\n\n";
		print FPR "[Input]\n";
		print FPR ">sample1\n";
}

#---------------------------------------------------------------------
###### Printing the last part of the parameter file of TAMA
sub PRINT_PARAM2{
		print FPR "[Preprocessing]\n";
		print FPR "\$TRIMMOMATIC-RUN=true\n";
		print FPR "\$TRIMMOMATIC-OPTION=\n";
		print FPR "\$BAYESHAMMER-RUN=true\n\n";
}
#---------------------------------------------------------------------

sub HELP{
		my $src = basename($0);
		print "\nUsage: $src [options] -1 <pe1> -2 <pe2> (or -s <se>) -o <out_dir>\n";
		print "\n== MANDATORY \n";
		print "-s\t<se>\tFile with unpaired reads [incompatible with -1 and -2]\n";
		print "-1\t<pe1>\tFile with #1 mates (paired 1) [incompatible with -s]\n";
		print "-2\t<pe2>\tFile with #2 mates (paired 2) [incompatible with -s]\n";
		print "-p|-core\t<integer>\tthe number of threads (default: 1)\n";
		print "-o\t<out_dir>\tOutput directory (default: ./ref_dir)\n";
		exit();
}
