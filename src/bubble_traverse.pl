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
		my @mapped_r = split(/,/,$reads);
## Counting 
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

my %hs_matrix = ();
my %hs_informative_read = ();
my %hs_bubble_count = ();
open(FBG, "$in_bubble_g");
while(<FBG>){
		chomp;
		if($_ eq ""){
				next;
		}
		my ($b_id, $container, $pos_z, $pos, $base, $r_count, $reads) = split(/\t/,$_);
		my @mapped_r = split(/,/,$reads);
		my @informative = ();
		my @useless = ();
		my ($bb, $bubble_num, $sub_id) = split(/\W/,$b_id);
		$hs_bubble_count{$container} = $bubble_num;
		for(my $i = 0; $i <= $#mapped_r; $i++){
				if($hs_read{$container}{$mapped_r[$i]} == 1){
						push(@useless, $mapped_r[$i]);
				}
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

print STDERR "#Container\tReadID\tlink_start_pos\tlink_end_pos\tfirst_continuous_link_size\tall_link_size\tBubble_link\n";
foreach my $container (sort keys(%hs_read)){
		my %hs_bubble_combi = ();
		my %hs_bubble_end = (); # added at Mar 16, 2021
		my @mapped_read = sort {$a<=>$b} keys(%{$hs_informative_read{$container}});
		foreach my $read (@mapped_read){
				my $link_start = 0;
				my $start_flag = 0;
				my $length_first_continuos_link = 0;
				my $length_flag = 0;
				my $length_link = 0;
				my @str_link = ();
				$hs_bubble_end{$read} = $hs_bubble_count{$container};
				for(my $i = 1; $i <= $hs_bubble_count{$container}; $i++){
						if(!exists($hs_matrix{$container}{$i}{$read})){
								push(@str_link,0);
								if($length_flag == 1){
										$length_flag = 2;
								}
						}
						else{
								$length_link += 1;
								$hs_bubble_end{$read} = $i;
								push(@str_link,$hs_matrix{$container}{$i}{$read});
								if($start_flag == 0 && $hs_matrix{$container}{$i}{$read} != 0){
										$link_start = $i;
										$start_flag = 1;
										$length_flag = 1;
								}
								if($length_flag == 2){ next; }
								else{
										$length_first_continuos_link += 1;
								}
						}
				}
				my $str_link = join(":", @str_link);
				print STDERR "$container\t$read\t$link_start\t$hs_bubble_end{$read}\t$length_first_continuos_link\t$length_link\t$str_link\n";
				$hs_bubble_combi{$link_start}{$length_first_continuos_link}{$length_link}{$read} = $str_link;
		}
#		print STDERR "\n";
######################################################### Find consensus
		my %hs_consensus = ();
		my $consensus_id = 0;
		foreach my $l_start (sort {$a<=>$b} keys(%hs_bubble_combi)){
				foreach my $f_link_len (sort {$b<=>$a} keys(%{$hs_bubble_combi{$l_start}})){
						foreach my $all_link_len (sort {$b<=>$a} keys(%{$hs_bubble_combi{$l_start}{$f_link_len}})){
								foreach my $r (sort {$a<=>$b} keys(%{$hs_bubble_combi{$l_start}{$f_link_len}{$all_link_len}})){
										my $cur_combi = $hs_bubble_combi{$l_start}{$f_link_len}{$all_link_len}{$r};
										if($consensus_id == 0){
												$consensus_id = 1;
												$hs_consensus{$consensus_id} = $cur_combi;
												print STDERR "#LOOP:$container\t$r\t-\t-\t$cur_combi\tinitial-$consensus_id\n";
										}
										else{
												my $is_update = 0;
												foreach my $cons_id (sort {$a<=>$b} keys(%hs_consensus)){
														my $is_equal = 1;
														my @cons = split(/:/,$hs_consensus{$cons_id});
														my @target = split(/:/, $cur_combi);
														my ($loop_start, $loop_end) = (($l_start-1), ($hs_bubble_end{$r}-1));
		#print "$hs_consensus{$cons_id}\n$cur_combi\n";
														my $overlaps = 0; ## Added at 4/22 02:27
														for(my $p = $loop_start; $p <= $loop_end; $p++){
																if($cons[$p] == $target[$p] && $cons[$p] != 0){ $overlaps += 1; } ## Added at 4/22 02:27
																if($cons[$p] != 0 && $target[$p] != 0 && $cons[$p] != $target[$p]){
																		$is_equal = -1;
																		last;
																}
														}
														if($is_equal == 1){ ## Consensus update
																if($overlaps == 0){ next; } ## Added at 4/22 02:27
																for(my $p = $loop_start; $p <= $loop_end; $p++){
																		if($cons[$p] == 0 && $target[$p] != 0){
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
												if($is_update == 0){
														$consensus_id += 1;
														$hs_consensus{$consensus_id} = $cur_combi;
														print STDERR "#LOOP:$container\t$r\t-\t-\t$cur_combi\tinitial-$consensus_id\n";
												}
										}
								}
						}
				}
		}
		foreach my $cons_id (sort {$a<=>$b} keys(%hs_consensus)){
				print "$container\t$cons_id\t$hs_consensus{$cons_id}\n";
		}
#		print "\n";
}
