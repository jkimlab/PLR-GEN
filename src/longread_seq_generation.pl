#!/usr/bin/perl
use strict;
use warnings;

#my $container_list = $ARGV[0]; # raw/merged.l1.c1.bed
my $in_normal_node = $ARGV[0]; # NORMAL merged
my $in_bubble_node = $ARGV[1]; # BUBBLE.graph.out
my $in_bubble_path = $ARGV[2]; # consensus_bubble_path.txt
my $min_length = $ARGV[3]; # default : 100
my $out_dir = $ARGV[-1];

`mkdir -p $out_dir`;
my %hs_node_seq = ();
my %hs_bubble_path = ();
my %hs_container_list = ();

################################# Normal nodes
open(FNNODE, "$in_normal_node");
while(<FNNODE>){
		chomp;
		my ($container, $start, $end, $seq) = split(/\t/,$_);
		$hs_node_seq{$container}{$start}{N} = $seq;
		if(!exists($hs_container_list{$container})){
				$hs_container_list{$container} = 1;
		}
}
close(FNNODE);
################################# Bubble nodes
open(FBNODE, "$in_bubble_node");
while(<FBNODE>){
		chomp;
		if($_ eq ""){ next; }
		my ($bubble_id, $container, $start, $end, $seq, $count, $reads) = split(/\t/,$_);
		$bubble_id = substr($bubble_id,7);
		$hs_node_seq{$container}{$start}{$bubble_id} = $seq;
		if(!exists($hs_container_list{$container})){
				$hs_container_list{$container} = 1;
		}
}
close(FBNODE);
################################# Bubble paths
open(FBPATH, "$in_bubble_path");
while(<FBPATH>){
		chomp;
		my ($container, $id, $path) = split(/\t/,$_);
		$hs_bubble_path{$container}{$id} = $path;
}
close(FBPATH);
open(FOUT, ">$out_dir/pseudo-long_read.fa");
################################# Sequence generation with bubble path
print STDERR "# With consensus bubble path\n";
my @container_with_bubble_path = (sort keys(%hs_bubble_path));
foreach my $container (@container_with_bubble_path){
		print STDERR ">$container\n";
		foreach my $id (sort {$a<=>$b} keys(%{$hs_bubble_path{$container}})){
				my @path = split(/:/,$hs_bubble_path{$container}{$id});
				print STDERR "Path $id : @path\n";
				my $cur_plr_seq = "";
				my $bubble_num = 0;
				my $p = -1;
				foreach my $start (sort {$a<=>$b} keys(%{$hs_node_seq{$container}})){
						my @node_type = keys(%{$hs_node_seq{$container}{$start}});
						if(scalar(@node_type) == 1){ ## normal node
								$cur_plr_seq .= $hs_node_seq{$container}{$start}{$node_type[0]};
						}
						else{
								$p++;
								my $selected_bubble = $path[$p];
								$bubble_num += 1;
								print STDERR "$bubble_num-$selected_bubble ";
								if($selected_bubble == 0){
										$cur_plr_seq .= "N";
										print STDERR "(N) ";
								}
								else{
										my $cur_bubble_id = "$bubble_num-$selected_bubble";
										$cur_plr_seq .= $hs_node_seq{$container}{$start}{$cur_bubble_id};
										print STDERR "($hs_node_seq{$container}{$start}{$cur_bubble_id}) ";
								}
						}
				}
				print STDERR "\n";
				if(length($cur_plr_seq) < $min_length){
						print STDERR "## $container:$id -> skipped (too short)\n";
				}
				else{
						print FOUT ">$container:$id\n$cur_plr_seq\n";
				}
		}
		delete($hs_container_list{$container});
}
################################# Sequence generation no bubble path
print STDERR "# No consensus bubble path\n";
foreach my $container (sort keys(%hs_container_list)){
		print STDERR ">$container\n";
		my $cur_plr_seq = "";
		foreach my $start (sort {$a<=>$b} keys(%{$hs_node_seq{$container}})){
				my @node_type = keys(%{$hs_node_seq{$container}{$start}});
				if(scalar(@node_type) == 1){ ## normal node
						$cur_plr_seq .= $hs_node_seq{$container}{$start}{$node_type[0]};
				}
				else{
						$cur_plr_seq .= "N";
				}
		}
		if(length($cur_plr_seq) < $min_length){
				print STDERR "## $container:1 -> skipped (too short)\n";
		}
		else{
				print FOUT ">$container:1\n$cur_plr_seq\n";
		}
}
