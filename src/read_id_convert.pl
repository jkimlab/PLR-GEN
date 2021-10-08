#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';

my $in_read1 = $ARGV[0];
my $in_read2 = $ARGV[1];
my $out_dir = $ARGV[2]; ## run_plr of $out_dir/tmp

if($in_read1 =~ /.gz$/){ open(FREAD1, "gunzip -c $in_read1 |"); }
else{ open(FREAD1, "$in_read1"); }
open(FIDMAP, ">$out_dir/converted_id.map");
open(FOUT1, ">$out_dir/read_1.fq");

my @read_ids = ();
my $num_of_read1 = 0;
my $num_of_read2 = 0;
my $line = -1;
my $read_num = 0;
while(<FREAD1>){
		chomp;
		$line++;
		if($line%4 == 0){
				if($_ !~ /^@/){
						print STDERR "There is missing line in the $in_read1\n";
						print "Wrong";
						exit;
				}
				$read_num++;
				$num_of_read1++;
				my $cur_read_id = "$read_num/1";
				my $original_id = substr($_, 1);
				print FIDMAP "$original_id\t$cur_read_id\n";
				print FOUT1 "\@$cur_read_id\n";
				push(@read_ids, $read_num);
		}
		else{
				print FOUT1 "$_\n";
		}
}
close(FREAD1);
close(FOUT1);

if($in_read2 =~ /.gz$/){ open(FREAD2, "gunzip -c $in_read2 |"); }
else{ open(FREAD2, "$in_read2"); }
open(FOUT2, ">$out_dir/read_2.fq");
$line = -1;
while(<FREAD2>){
		chomp;
		$line++;
		$num_of_read2++;
		if($line%4 == 0){
				if($_ !~ /^@/){
						print STDERR "There is missing line in the $in_read2\n";
						print "Wrong";
						exit;
				}
				my $cur_read_id = shift(@read_ids);
				my $original_id = substr($_, 1);
				print FIDMAP "$original_id\t$cur_read_id/2\n";
				print FOUT2 "\@$cur_read_id/2\n";
		}
		else{
				print FOUT2 "$_\n";
		}
}
close(FREAD2);
close(FOUT2);
close(FIDMAP);
if($#read_ids != -1){ 
		print "Wrong\n";
}


