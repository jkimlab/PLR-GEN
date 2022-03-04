#!/usr/bin/env perl
use strict;
use warnings;
use Scalar::Util qw(looks_like_number);
use FindBin qw($Bin);

my $in_mpileup = $ARGV[0];
my $min_length = $ARGV[1];
my $min_cnt = $ARGV[2];
my $out_dir = $ARGV[-1];

#---------------------------------------------------------------------
###### Preprocessing 1 : extracting aligned columns with alignment depth >= min_cutoff
$min_length = $min_length/10;
`mkdir -p $out_dir`;
`awk '\$4>=$min_cnt {print \$1 \"\t\" \$2-1 \"\t\" \$2 \"\t\" \$3 \"\t\" \$4 \"\t\" \$5 \"\t\" \$6 \"\t\" \$7}' $in_mpileup > $out_dir/mapped.bed`;
if(-z "$out_dir/mapped.bed"){
		`touch $out_dir/.no_mapping`;
		exit();
}
#---------------------------------------------------------------------
###### Preprocessing 2 : creating PLR containers
`cut -f1,2,3 $out_dir/mapped.bed > $out_dir/pos.bed`; 
`$Bin/../bin/bedtools sort -i $out_dir/pos.bed > $out_dir/sort.pos.bed`;
`$Bin/../bin/bedtools merge -i $out_dir/sort.pos.bed | awk '\$3-\$2>=$min_length' - > $out_dir/container_list.bed`;
`$Bin/../bin/bedtools intersect -wo -a $out_dir/container_list.bed -b $out_dir/mapped.bed > $out_dir/mapping_info_in_container.txt`;

#---------------------------------------------------------------------
###### Making PLR graphs
open(FOGRAPH, ">$out_dir/graph.out");
open(FCONTAINER, "$out_dir/mapping_info_in_container.txt");
my %hs_list = ();
my $saved_node = "";
my $pre_container = "";
my $coord = 0;
my $bubble_id = 0;

######### converting alignment columns to nodes of PLR graphs
while(<FCONTAINER>){ # for each base position
		chomp;
## parsing information of an alignment 
		$coord += 1;
		my ($c1, $c2, $c3, $m1, $t1, $pos, $r_base, $count, $aln, $mapq, $reads, $t2) = split(/\t/,$_);
		my $container_name = "$c1:$c2-$c3";
		if($pre_container ne $container_name){
				$pre_container = $container_name;
				$coord = 1;
				$bubble_id=0;
		}
		my $processed_aln = $aln;
		$processed_aln =~ tr/,acgt/.ACGT/;
		my @aln = split(//,$processed_aln);
		my @g_node = ();
		my ($cur_base, $next_base) = ("", "");
		for(my $i = 0; $i <= $#aln; $i++){

## parsing aligned bases
				if($aln[$i] eq "^"){
						$i += 1; #skip next character (mapping quality)
						next;
				}
				if($aln[$i] eq '$'){
#						next;
				}
				elsif($aln[$i] eq "."){ $cur_base = $r_base; }
				elsif($aln[$i] eq "A" || $aln[$i] eq "C" || $aln[$i] eq "G" || $aln[$i] eq "T"){ $cur_base = $aln[$i]; }
				else{ next; }
				if($i == $#aln){ $next_base = ""; }
				else{ $next_base = $aln[$i+1]; }

				if($next_base eq "."){ push(@g_node, $cur_base); next; }
				elsif($next_base eq "A" || $next_base eq "C" || $next_base eq "G" || $next_base eq "T" || $next_base eq ""){ push(@g_node, $cur_base); next; }
				elsif($next_base eq "+" || $next_base eq "-"){
						my $cur_node = $cur_base;
						my $insert_size = "";
						my $extracting_start = 0;
						for(my $n = $i+2; $n <= $#aln; $n++){
								if(looks_like_number($aln[$n])){ $insert_size .= $aln[$n]; }
								else{
										$extracting_start = $n;
										$insert_size = int($insert_size);
										last;
								}
						}
						$i = $extracting_start+$insert_size-1;
						my $insert_base = substr($processed_aln, $extracting_start, $insert_size);
						if($next_base eq "+"){
								push(@g_node, "$cur_base$insert_base");
						}
						elsif($next_base eq "-" ){
								push(@g_node, "$cur_base");
						}
				}
		}
		my $depth = scalar(@g_node);
		if($depth < $min_cnt){ next; }

## counting mapping depth for each base (A, C, G, and T)
		my %hs_cnt = ();
		my %hs_read = ();
		my @read = split(/,/,$reads);
		for(my $i = 0; $i <= $#g_node; $i++){
				my $node = $g_node[$i];
				if(!exists($hs_cnt{$node})){
						$hs_cnt{$node} = 1; 
						$hs_read{$node} = $read[$i];
				}
				else{
						$hs_cnt{$node} += 1;
						$hs_read{$node} .= ",$read[$i]";
				}
		}

## creating normal and bubble nodes
		my @count = sort {$b<=>$a} values(%hs_cnt);
		my $s_coord = $coord-1;
		if($count[0] == scalar(@g_node)){
				print FOGRAPH "NORMAL\t$container_name\t$s_coord\t$coord\t$g_node[0]\n";
		}
		else{
				foreach my $b_node (keys(%hs_cnt)){
						my @reads = split(/,/,$hs_read{$b_node});
						@reads = sort {$a<=>$b} @reads;
						my $mapped_read = join(",",@reads);
						$bubble_id += 1;
						print FOGRAPH "BUBBLE\t$container_name\t$s_coord\t$coord\t$b_node\t$hs_cnt{$b_node}\t$mapped_read\n";
				}
		}
}
close(FCONTAINER);

#---------------------------------------------------------------------

