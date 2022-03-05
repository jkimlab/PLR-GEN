#!/usr/bin/env perl
use strict;
use warnings;

my $in_bubble_g = shift;

#---------------------------------------------------------------------
###### Storing information of mapped reads
my %hs_read = ();
open(FBG, "$in_bubble_g");
while(<FBG>){
		chomp;
		if($_ eq ""){ next; }
		my ($b_id, $container, $pos_z, $pos, $base, $r_count, $reads) = split(/\t/,$_);
		my @mapped_r = split(/,/,$reads); # all mapped read IDs are stored in array @mapped_r
## Counting the number of overlapped bubbles in a read
		for(my $i = 0; $i <= $#mapped_r; $i++){
				if(!exists($hs_read{$container}{$mapped_r[$i]})){
						$hs_read{$container}{$mapped_r[$i]} = 1;
				}
				else{
						$hs_read{$container}{$mapped_r[$i]} += 1;
				}
		}
}
close(FBG);
#---------------------------------------------------------------------
###### Finding and storing informative reads
my %hs_matrix = ();
my %hs_informative_read = ();
my %hs_bubble_count = ();
open(FBG, "$in_bubble_g");
while(<FBG>){
		chomp;
		if($_ eq ""){ next; }
		my ($b_id, $container, $pos_z, $pos, $base, $r_count, $reads) = split(/\t/,$_);
		my @mapped_r = split(/,/,$reads);
		my @informative = ();
		my @useless = ();
		my ($bb, $bubble_num, $sub_id) = split(/\W/,$b_id);
		$hs_bubble_count{$container} = $bubble_num;
## Filtering out reads that can not provide linking information of bubbles 
#### Reads overlapped only one bubble will be ignored 
		for(my $i = 0; $i <= $#mapped_r; $i++){ 
				if($hs_read{$container}{$mapped_r[$i]} == 1){
						push(@useless, $mapped_r[$i]);
				}
#### Reads overlapped two or more bubbles will be stored
				else{
						push(@informative, $mapped_r[$i]);
						$hs_matrix{$container}{$bubble_num}{$mapped_r[$i]} = $sub_id;
						$hs_informative_read{$container}{$mapped_r[$i]} = 1;
				}
		}
		my $infor = join(",",@informative);
		my $nouse = join(",",@useless);
}
close(FBG);

#---------------------------------------------------------------------
###### Finding bubble combinations from informative reads
print STDERR "#Container\tReadID\tlink_start_pos\tlink_end_pos\tfirst_continuous_link_size\tall_link_size\tBubble_link\n";
## For each container,
foreach my $container (sort keys(%hs_read)){
		my %hs_bubble_combi = ();
		my %hs_bubble_end = (); 
		my @mapped_read = sort {$a<=>$b} keys(%{$hs_informative_read{$container}});
## For each informative read,
		foreach my $read (@mapped_read){
				my $link_start = 0;
				my $start_flag = 0;
				my $length_first_continuos_link = 0;
				my $length_flag = 0;
				my $length_link = 0;
				my @str_link = ();
				$hs_bubble_end{$read} = $hs_bubble_count{$container};
## Searching which bubbles are overlapped in a read
#### For ith bubble,
				for(my $i = 1; $i <= $hs_bubble_count{$container}; $i++){
						if(!exists($hs_matrix{$container}{$i}{$read})){ ## ith bubble is not overlapped
								push(@str_link,0);
								if($length_flag == 1){
										$length_flag = 2;
								}
						}
						else{ ## ith bubble is overlapped
								$length_link += 1;
								$hs_bubble_end{$read} = $i;
								push(@str_link,$hs_matrix{$container}{$i}{$read});
								if($start_flag == 0 && $hs_matrix{$container}{$i}{$read} != 0){ ## check if it is the first overlapped bubble of this read
										$link_start = $i;
										$start_flag = 1;
										$length_flag = 1;
								}
								if($length_flag == 2){ next; }  ## check if this is the first overlapped bubble block (a group of contiguous bubbles) of this read
								else{
										$length_first_continuos_link += 1;
								}
						}
				}
				my $str_link = join(":", @str_link);  ## $str_link is the final combination of bubbles supported by this read 
				print STDERR "$container\t$read\t$link_start\t$hs_bubble_end{$read}\t$length_first_continuos_link\t$length_link\t$str_link\n";
				$hs_bubble_combi{$link_start}{$length_first_continuos_link}{$length_link}{$read} = $str_link;
		}
		
######################################################### Find consensus
		my %hs_consensus = ();
		my $consensus_id = 0;
		foreach my $l_start (sort {$a<=>$b} keys(%hs_bubble_combi)){  # (1) Start with the first overlapped bubble 
				foreach my $f_link_len (sort {$b<=>$a} keys(%{$hs_bubble_combi{$l_start}})){  # (2) among them, start with the longest overlapped bubble block
						foreach my $all_link_len (sort {$b<=>$a} keys(%{$hs_bubble_combi{$l_start}{$f_link_len}})){ # (3) among them, start with the combination with the most overlapped bubbles
								foreach my $r (sort {$a<=>$b} keys(%{$hs_bubble_combi{$l_start}{$f_link_len}{$all_link_len}})){
										my $cur_combi = $hs_bubble_combi{$l_start}{$f_link_len}{$all_link_len}{$r};
										if($consensus_id == 0){  ## Initialization
												$consensus_id = 1;
												$hs_consensus{$consensus_id} = $cur_combi;
												print STDERR "#LOOP:$container\t$r\t-\t-\t$cur_combi\tinitial-$consensus_id\n";
										}
										else{
												my $is_update = 0;
												foreach my $cons_id (sort {$a<=>$b} keys(%hs_consensus)){
														my $is_equal = 1;
														my @cons = split(/:/,$hs_consensus{$cons_id});  ## the combination of this consensus of combination of bubbles
														my @target = split(/:/, $cur_combi); ## the combination of bubble of this read
														my ($loop_start, $loop_end) = (($l_start-1), ($hs_bubble_end{$r}-1)); 
														my $overlaps = 0;
														for(my $p = $loop_start; $p <= $loop_end; $p++){
																if($cons[$p] == $target[$p] && $cons[$p] != 0){ $overlaps += 1; } ## linked bubbles are the same at this position
																if($cons[$p] != 0 && $target[$p] != 0 && $cons[$p] != $target[$p]){ ## linked bubbles are not the same at this position (-> break this loop)
																		$is_equal = -1;
																		last;
																}
														}
														if($is_equal == 1){ ## Consensus update
																if($overlaps == 0){ next; } ## there is no overlapped position with the information of linked bubbles
																for(my $p = $loop_start; $p <= $loop_end; $p++){
																		if($cons[$p] == 0 && $target[$p] != 0){ ## fill the gap of consensus with the linked bubble of this combination
																				$cons[$p] = $target[$p];
																		}
																}
																my $updated_cons = join(":",@cons);
																$hs_consensus{$cons_id} = $updated_cons;
																$is_update = 1;
																print STDERR "#LOOP:$container\t$r\t$loop_start\t$loop_end\t$cur_combi\tupdated-$cons_id\n";
																last;
														}
												}
												if($is_update == 0){  ## make a new consensus 
														$consensus_id += 1;
														$hs_consensus{$consensus_id} = $cur_combi;
														print STDERR "#LOOP:$container\t$r\t-\t-\t$cur_combi\tinitial-$consensus_id\n";
												}
										}
								}
						}
				}
		}
#---------------------------------------------------------------------
###### Finishing for this conatiner
		foreach my $cons_id (sort {$a<=>$b} keys(%hs_consensus)){
				print "$container\t$cons_id\t$hs_consensus{$cons_id}\n";
		}
}
