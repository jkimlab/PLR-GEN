#!/usr/bin/perl
use strict;
use warnings;

my $in_graph = shift;

my %hs_stat = ();
open(FGR, "$in_graph");
while(<FGR>){
		chomp;
		my @t = split(/\t/,$_);
		if($t[0] !~ /^BUBBLE/){ next; }
		if(!exists($hs_stat{$t[1]}{$t[3]})){
				$hs_stat{$t[1]}{$t[3]} = "$t[5]";
		}
		else{
				$hs_stat{$t[1]}{$t[3]} .= ",$t[5]";
		}
}
close(FGR);
foreach my $container (sort keys(%hs_stat)){
		my $total_bubbles = 0;
		foreach my $pos (sort {$a<=>$b} keys(%{$hs_stat{$container}})){
				my @counts = split(",",$hs_stat{$container}{$pos});
				@counts = sort {$b<=>$a} @counts;
				my $n_counts = scalar(@counts);
				my $total_counts = 0;
				for(my $i = 0; $i <= $#counts; $i++){
						$total_counts += $counts[$i];
				}
				print "$container\t$pos\t$n_counts\t$total_counts";
				for(my $i = 0; $i <= $#counts; $i++){
						#my $prop = sprintf("%.3f",($counts[$i]/$total_counts));
						my $prop = ($counts[$i]/$total_counts);
						print "\t$prop";
#						print "\t$counts[$i]";
				}
				print "\n";
				$total_bubbles += 1;
		}
		print STDERR "$container\t$total_bubbles\n";
}
