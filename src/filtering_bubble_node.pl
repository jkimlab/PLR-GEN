#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use Math::Round;

## Arguments
my $in_graph = $ARGV[0];
my $in_bubble_stat = $ARGV[1];
my $in_cutoff_perc = $ARGV[2];
my $in_cutoff_minor = $ARGV[3];
my $out_dir = $ARGV[-1];

`mkdir -p $out_dir`;

my @arr_read_count = ();
my @arr_minor_freq = ();
my %hs_filtout = ();
my %hs_BtoN = ();
my %hs_n_node = ();
#---------------------------------------------------------------------
###### Storing information of bubbles
if(!-z "$in_bubble_stat"){ ## If there is no bubble, this process will be skipped
		open(FB, "$in_bubble_stat");
		while(<FB>){
				chomp;
				my ($container, $pos, $bubble_count, $read_count, $max_freq, @minor_freq) = split(/\t/,$_);
				push(@arr_read_count, $read_count);
		}
		close(FB);
		
#---------------------------------------------------------------------
###### Converting bubbles to normal nodes 
		@arr_read_count = sort {$a<=>$b} @arr_read_count;
		
## Filtering 1: bubbles with low read depth
#### -> bubbles with low read depth < cutoff will be converted to normal node with "N" base
		my $n_bubble = scalar(@arr_read_count);
		my $b_pi = int($n_bubble*($in_cutoff_perc/100))-1;
		my $cutoff_n = $arr_read_count[$b_pi];
		my $cutoff_f = (1/$cutoff_n); ## -> the second filtering cutoff
		print STDERR "# read depth cutoff = $cutoff_n\n";
		print STDERR "# filtered out bubbles\n";

## Filtering 2: pruning bubble nodes with relatively low read depth in a bubble
		open(FSTAT, "$in_bubble_stat");
		while(<FSTAT>){
				chomp;
				my @t = split(/\t/,$_);
				my $survive = 1;
				my $bubble_to_normal = 1;
				my $count = $t[3];
				my $max = $t[4];
				if($count < $cutoff_n){ $bubble_to_normal = 0; }
				else{
						for(my $i = 5; $i <= $#t; $i++){
								my $prop = ($t[$i]/$max);
								if($in_cutoff_minor == 1 && $prop == 1){
										$survive = 1;
										$bubble_to_normal = 0; 
										last;
								}
								if($prop >= $in_cutoff_minor){ $survive += 1; }
						}
				}
				if($survive == 1){ # when bubble count is too low, $survive must be 1
						my $m_count = round($max*$t[3]); 
						if($bubble_to_normal == 1){ # This bubble will be changed to normal node with its base of the major bubble node
								$hs_BtoN{$t[0]}{$t[1]}{$m_count} = 1;
								print STDERR "BTON-1\t$t[0]\t$t[1]\t$m_count\n";
						}
						elsif($bubble_to_normal == 0){ # This bubble has to be changed to normal node with base "N"
								$hs_BtoN{$t[0]}{$t[1]}{$m_count} = 0;
								print STDERR "BTON-0\t$t[0]\t$t[1]\t$m_count\n";
						}
						for(my $i = 5; $i <= $#t; $i++){
								my $count = round($t[3]*$t[$i]);
								$hs_filtout{$t[0]}{$t[1]}{$count} = 1;
								print STDERR "FILT-O$t[0]\t$t[1]\t$count\n";
						}
				}
				else{
						for(my $i = 5; $i <= $#t; $i++){
								my $prop = ($t[$i]/$max);
								if($prop < $in_cutoff_minor){
										my $count = round($t[3]*$t[$i]);
										$hs_filtout{$t[0]}{$t[1]}{$count} = 1;
										print STDERR "FILT-O\t$t[0]\t$t[1]\t$count\n";
								}
						}
				}
		}
		close(FSTAT);
}
#---------------------------------------------------------------------
###### Making final graph
my %hs_bubble = ();
my %hs_sub_id = ();

open(FNOUT, ">$out_dir/NORMAL.graph.out");
open(FBOUT, ">$out_dir/BUBBLE.graph.out");
open(FGRAPH, "$in_graph");
my $b_id = 0;
while(<FGRAPH>){
		chomp;
		my @t = split(/\t/,$_);
		if($_ =~ /^NORMAL/){
				print FNOUT "$t[1]\t$t[2]\t$t[3]\t$t[4]\n";
				print "$t[1]\t$t[2]\t$t[3]\t$t[4]\n";
		}
		else{
				if(exists($hs_BtoN{$t[1]}{$t[3]}{$t[5]})){
						if($hs_BtoN{$t[1]}{$t[3]}{$t[5]} == 1){
								print FNOUT "$t[1]\t$t[2]\t$t[3]\t$t[4]\n";
								print "$t[1]\t$t[2]\t$t[3]\t$t[4]\n";
						}
						else{
								print FNOUT "$t[1]\t$t[2]\t$t[3]\tN\n";
								print "$t[1]\t$t[2]\t$t[3]\tN\n";
						}
				}
				else{
						if(!exists($hs_filtout{$t[1]}{$t[3]}{$t[5]})){
								if(!exists($hs_bubble{$t[1]})){
										print FBOUT "\n";
										$hs_bubble{$t[1]} = 1;
										$hs_sub_id{$t[1]}{$t[2]} = 1;
								}
								else{
										if(!exists($hs_sub_id{$t[1]}{$t[2]})){
												print FBOUT "\n";
												$hs_bubble{$t[1]} += 1; 
												$hs_sub_id{$t[1]}{$t[2]} = 1;
										}
										else{
												$hs_sub_id{$t[1]}{$t[2]} += 1;
										}
								}
										
								print FBOUT "BUBBLE:$hs_bubble{$t[1]}-$hs_sub_id{$t[1]}{$t[2]}\t$t[1]\t$t[2]\t$t[3]\t$t[4]\t$t[5]\t$t[6]\n";
								print "BUBBLE:$hs_bubble{$t[1]}-$hs_sub_id{$t[1]}{$t[2]}\t$t[1]\t$t[2]\t$t[3]\t$t[4]\t$t[5]\t$t[6]\n";
						}
				}
		}
}
close(FGRAPH);
