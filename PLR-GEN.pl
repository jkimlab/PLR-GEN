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
my $out_dir = "./PR_out";
my $in_cpu = 1;
my $ref_sampling;
my $mapping_quality = 20;
my $bubble_depth_cutoff = 1;
my $min_plr_length = 100;
my $min_aln_count = 1;
my $minor_bubble_cutoff = 0.5;
my $run_tama;
my $temp_save;
my $help;


## parameters
my $options = GetOptions(
		"s=s" => \$read_s, # mandatory (incompatible with -1 and -2)
		"1=s" => \$read1, # mandatory (incompatible with -s)
		"2=s" => \$read2, # mandatory (incompatible with -s)
		"r|ref=s" => \$in_ref_list, # mandatory (incompatible with -tama)
		"tama" => \$run_tama, # mandatory (incompatible with -ref)
		"sampling=f" => \$ref_sampling, # default: off;
		"o=s" => \$out_dir, # default: ./PR.out
		"p|core=i" => \$in_cpu, # default: 1
		"q|mapq=i" => \$mapping_quality, # default: 20
		"d|min_depth=i" => \$bubble_depth_cutoff, # default: 1
		"b|minor_bubble=f" => \$minor_bubble_cutoff, # default: 0.5
		"l|min_length=i" => \$min_plr_length, # default: 100
		"c|min_count=i" => \$min_aln_count,
		"t|temp" => \$temp_save, # all temporary file will be saved
		"h|help" => \$help # print help
);

##### Parameter check
if(!defined($options)){ HELP(); }
if(defined($help)){ HELP(); }
if(!defined($in_ref_list) && !defined($run_tama)){ HELP(); }
if(defined($in_ref_list) && defined($run_tama)){ HELP(); }
if(defined($read_s) && defined($read1)){
		HELP();
}
if(defined($read1) && !defined($read2)){
		HELP();
}
if(defined($read2) && !defined($read1)){
		HELP();
}
`mkdir -p $out_dir/data/ref $out_dir/raw_plrs`;
$out_dir = abs_path("$out_dir");
open(FLOG, ">$out_dir/log.PSEUDO-LONG_READ_GEN.txt");
print STDERR "= Used parameters\n";
print FLOG "= Used parameters\n";
if(defined($read_s)){
		print STDERR "== READ : $read_s\n";
}
else{
		print STDERR "== READ1 : $read1\n";
		print FLOG "== READ1 : $read1\n";
		print STDERR "== READ2 : $read2\n";
		print FLOG "== READ2 : $read2\n";
}
print STDERR "== REF LIST : $in_ref_list\n";
print FLOG "== REF LIST : $in_ref_list\n";
print STDERR "== CORES : $in_cpu\n";
print FLOG "== CORES : $in_cpu\n";
print STDERR "== MAPPING QUALITY : $mapping_quality\n";
print FLOG "== MAPPING QUALITY : $mapping_quality\n";
print STDERR "== MIN PSEUDO-LONG READ LENGTH : $min_plr_length\n";
print FLOG "== MIN PSEUDO-LONG READ LENGTH : $min_plr_length\n";
print STDERR "== MIN MAPPED READ DEPTH : $min_aln_count\n";
print FLOG "== MIN MAPPED READ DEPTH : $min_aln_count\n";
print STDERR "== MIN READ DEPTH OF BUBBLE : $bubble_depth_cutoff\n";
print FLOG "== MIN READ DEPTH OF BUBBLE : $bubble_depth_cutoff\n";
if(defined($temp_save)){
		print STDERR "== all intermediated files will be saved.\n";
		print FLOG "== all intermediated files will be saved.\n";
}
else{
		print STDERR "== all intermediated files will be removed.\n";
		print FLOG "== all intermediated files will be removed.\n";
}
print STDERR "== OUTPUT DIR : $out_dir\n\n";
print FLOG "== OUTPUT DIR : $out_dir\n\n";
#---------------------------------------------------------------------

##### Preprocessing
###### checking read pairs (if input is paired-end) converting read IDs
print STDERR "= Preprocessing\n";
print FLOG "= Preprocessing\n";
print STDERR "== READ ID CONVERTING\n";
print FLOG "== READ ID CONVERTING\n";
if(defined($read_s)){
		`$Bin/src/read_id_convert_single.pl $read_s $out_dir/data`;
		$read1 = "$out_dir/data/read.fq";
}
else{
		`$Bin/src/read_id_convert.pl $read1 $read2 $out_dir/data`;
		$read1 = "$out_dir/data/read_1.fq";
		$read2 = "$out_dir/data/read_2.fq";
}
#---------------------------------------------------------------------
###### preparing references
my %hs_ref = ();
my @arr_err = ();
my $check_all_ref = 0;
####### running tama (if -tama option is true)
if(defined($run_tama)){
		print STDERR "Running TAMA for preparing references\n";
		print FLOG "Running TAMA for preparing references\n";
		if(defined($read_s)){
				`$Bin/src/ref_preparation_TAMA.pl -s $read1 -p $in_cpu  -o $out_dir/data/tamaout`;
		}
		else{
				`$Bin/src/ref_preparation_TAMA.pl -1 $read1 -2 $read2 -p $in_cpu -o $out_dir/data/tamaout`;
		}
		$in_ref_list = "$out_dir/data/tamaout/ref_list.txt";
}
#---------------------------------------------------------------------

####### cheking reference files
print STDERR "== REFERENCE CHECK\n";
print FLOG "== REFERENCE CHECK\n";
open(FREF, "$in_ref_list");
while(<FREF>){
		chomp;
		my $ref_file = $_;
		if(!-f $ref_file){
				push(@arr_err, $ref_file);
		}
		else{
				$ref_file = abs_path("$ref_file");
				$hs_ref{$ref_file} = 1;
		}
}
close(FREF);
if(scalar(@arr_err) >= 1){
		print STDERR "Please check the reference file.\n";
		print FLOG "Please check the reference file.\n";
		foreach my $file (@arr_err){
				print STDERR "\t$file\n";
				print FLOG "\t$file\n";
		}
		exit();
}
#---------------------------------------------------------------------

my $num_of_all_refs = scalar(keys(%hs_ref));
print STDERR "Total number of references: $num_of_all_refs\n";
print FLOG "Total number of references: $num_of_all_refs\n";
my %hs_sample = ();
my @arr_all_ref_files = sort keys(%hs_ref);
####### random sampling for references (if -sampling option is used)
if(defined($ref_sampling)){
		print STDERR "Random sampling for references\n";
		print FLOG "Random sampling for references\n";
		my $num_sampling = int($num_of_all_refs*$ref_sampling);
		print STDERR "The number of used references: $num_sampling ($ref_sampling)\n";
		print FLOG "The number of used references: $num_sampling ($ref_sampling)\n";
		srand();
		my $i = 0;
		while(1){
				my $sample_index = int(rand($num_of_all_refs));
				if(!exists($hs_sample{$arr_all_ref_files[$sample_index]})){
						$hs_sample{$arr_all_ref_files[$sample_index]} = 1;
						print STDERR "REF $i : $arr_all_ref_files[$sample_index]\n";
						print FLOG "REF $i : $arr_all_ref_files[$sample_index]\n";
						$i++;
				}
				if($i > $num_sampling){ last; }
		}
		%hs_ref = %hs_sample;
		%hs_sample=();
}
#---------------------------------------------------------------------

##### Pseudo-long read generation
print STDERR "= Generating pseudo-long read sequences\n";
print FLOG "= Generating pseudo-long read sequences\n";
my @all_refs = sort keys(%hs_ref);
my @coms = ();
my $num_refs = scalar(@all_refs);
my $num_cpu = $in_cpu;
my $num_process = 1;
if($in_cpu >= 10){
		if($num_process > $num_refs){ $num_process = $num_refs; }
		else{ $num_process = 5; }
		$num_cpu = int($in_cpu/$num_process);
}

my $FM_lo = new Parallel::ForkManager($num_process);
for(my $i = 0; $i < scalar(@all_refs); $i++){
		print STDERR "ref $i : $all_refs[$i]\n";
		print FLOG "ref $i : $all_refs[$i]\n";
		if($all_refs[$i] =~ /.gz$/){ `gunzip -c $all_refs[$i] > $out_dir/data/ref/ref_$i.fa`; }
		else{
				`cp $all_refs[$i] $out_dir/data/ref/ref_$i.fa`;
		}
		my $lo_qid = $FM_lo -> start and next;
		DO_PLR_GEN("$out_dir/data/ref/ref_$i.fa","out$i",$i); # Generating PLR sequence using the ith reference
		$FM_lo -> finish;
}
$FM_lo -> wait_all_children;
#---------------------------------------------------------------------

##### Finishing
print STDERR "= Merging all pseudo-long reads\n";
print FLOG "= Merging all pseudo-long reads\n";
`cat $out_dir/raw_plrs/*.pseudo-long_read.fa > $out_dir/pseudo-long_read.fa`;
`gzip $out_dir/raw_plrs/*`;
print STDERR "= Finished\n";
print FLOG "= Finished\n";
#---------------------------------------------------------------------

sub HELP{
		my $src = basename($0);
		print "\nUsage: $src [options] -1 <pe1> -2 <pe2> (or -s <se>) -r <ref_list> -o <out_dir>\n";
		print "\n== MANDATORY \n";
		print "-s\t<se>\tFile with unpaired reads [incompatible with -1 and -2]\n";
		print "-1\t<pe1>\tFile with #1 mates (paired 1) [incompatible with -s]\n";
		print "-2\t<pe2>\tFile with #2 mates (paired 2) [incompatible with -s]\n";
		print "-r|-ref\t<ref_list>\tThe list of reference genome sequence files\n";
		print "-tama\tReference preparation using TAMA [incompatible with -r|-ref]\n";
		print "-sampling\t<proportion>\tproportion to random sampling for references (default: off, range: 0-1)\n";
		print "-o\t<out_dir>\tOutput directory (default: ./PR.out)\n";
		print "\n==Running and filtering options\n";
		print "-p|-core\t<integer>\tThe number of threads (default: 1)\n";
		print "-q|-mapq\t<integer>\tMinimum mapping quality (default: 20)\n";
		print "-l|-min_length\t<integer>\tCutoff of minimum length of pseudo-long reads (default: 100bp)\n";
		print "-c|-min_count\t<integer>\tCutoff of minimum mapping depth for each node (default: 1)\n";
		print "-d|-min_depth\t<integer>\ Cutoff of mapping depth of bubbles (default: 1, range: 0-100)\n";
		print "\t\t0: all bubbles are used.\n";
		print "\t\t1: bubbles with less than 1% mapping depth from mapping depth distribution of bubbles are converted to normal nodes.\n";
		print "\t\t100: all bubbles are converted to normal nodes\n";
		print "\n==Other options\n";
		print "-t|-temp\tIf you use -t option, all intermediate files are left.\n";
		print "\tPlease careful to use this option because it has to be needed very large space.\n";
		print "-h|-help\tPrint help page.\n";
		exit;
}

##### Main workflow of PLR-GEN
sub DO_PLR_GEN{
		my $ref_fa = shift(@_);
		my $plr_outdir = shift(@_);
		my $plr_id = shift(@_);
		print FLOG "== $plr_id ... started\n";
		`mkdir -p $out_dir/$plr_outdir/raw`;
###### mapping reads to reference
		`$Bin/bin/bowtie2-build -f --threads $num_cpu $ref_fa $out_dir/$plr_outdir/raw/index`;
		if(!defined($read2)){
				`$Bin/bin/bowtie2 -p $num_cpu -x $out_dir/$plr_outdir/raw/index -U $read1 | $Bin/bin/samtools sort -o $out_dir/$plr_outdir/raw/mapping.bam -T $out_dir/$plr_outdir/raw/tmp -`;
		}
		else{
				`$Bin/bin/bowtie2 -p $num_cpu -x $out_dir/$plr_outdir/raw/index -1 $read1 -2 $read2 | $Bin/bin/samtools sort -o $out_dir/$plr_outdir/raw/mapping.bam -T $out_dir/$plr_outdir/raw/tmp -`;
		}
#---------------------------------------------------------------------

####### PLR generation (PLR container construction and PLR sequence generation)
		`bash $Bin/src/runMakeGraph.sh $Bin/src $Bin/bin $out_dir/$plr_outdir/raw/mapping.bam $ref_fa $mapping_quality $bubble_depth_cutoff $min_plr_length $min_aln_count $minor_bubble_cutoff $out_dir/$plr_outdir $plr_id`;
#---------------------------------------------------------------------

####### Finishing
		if(-f "$out_dir/$plr_outdir/pseudo-long_read.fa"){
				`cp $out_dir/$plr_outdir/pseudo-long_read.fa $out_dir/raw_plrs/$plr_id.pseudo-long_read.fa`;
		}
		if(-f "$out_dir/$plr_outdir/pseudo-long_read.fa.gz"){
				`gunzip -c $out_dir/$plr_outdir/pseudo-long_read.fa.gz > $out_dir/raw_plrs/$plr_id.pseudo-long_read.fa`;
		}
		if(!defined($temp_save)){
				`rm -rf $out_dir/$plr_outdir`;
		}
		print FLOG "== $plr_id ... finished\n";
#---------------------------------------------------------------------
}
close(FLOG);
