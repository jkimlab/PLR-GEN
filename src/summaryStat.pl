#!/usr/bin/env perl
use strict;
use warnings;

my $in_bubble_g = $ARGV[0];

my @arr_read_count = ();
my @arr_minor_freq = ();
open(FB, "$in_bubble_g");
while(<FB>){
		chomp;
		my ($container, $pos, $bubble_count, $read_count, $max_freq, @minor_freq) = split(/\t/,$_);
		push(@arr_read_count, $read_count);
		push(@arr_minor_freq, @minor_freq);
}
close(FB);

@arr_read_count = sort {$a<=>$b} @arr_read_count;
@arr_minor_freq = sort {$a<=>$b} @arr_minor_freq;

my $n_bubble = scalar(@arr_read_count);
my $n_minor_freq = scalar(@arr_minor_freq);
print "percentile\tRead_count\tMinor_frequency\n";
print "min\t$arr_read_count[0]\t$arr_minor_freq[0]\n";
for(my $i = 1; $i <= 100; $i++){
		my $b_pi = int($n_bubble*($i/100))-1;
		my $f_pi = int($n_minor_freq*($i/100))-1;
		print "$i\t$arr_read_count[$b_pi]\t$arr_minor_freq[$f_pi]\n";
}
print "max\t$arr_read_count[-1]\t$arr_minor_freq[-1]\n";

=pod ## Filtering test output ! 
my $perc_5 = int($n_bubble*(5/100))-1;
my $n_cutoff = $arr_read_count[$perc_5];
print "=== cutoff : $n_cutoff\n";
@arr_read_count = ();
@arr_minor_freq = ();
open(FB, "$in_bubble_g");
while(<FB>){
		chomp;
		my ($container, $pos, $bubble_count, $read_count, $max_freq, @minor_freq) = split(/\t/,$_);
		if($read_count < $n_cutoff){ next; }
		push(@arr_read_count, $read_count);
		push(@arr_minor_freq, @minor_freq);
}
close(FB);

@arr_read_count = sort {$a<=>$b} @arr_read_count;
@arr_minor_freq = sort {$a<=>$b} @arr_minor_freq;

$n_bubble = scalar(@arr_read_count);
$n_minor_freq = scalar(@arr_minor_freq);
print "percentile\tRead_count\tMinor_frequency\n";
print "min\t$arr_read_count[0]\t$arr_minor_freq[0]\n";
for(my $i = 1; $i <= 100; $i++){
		my $b_pi = int($n_bubble*($i/100))-1;
		my $f_pi = int($n_minor_freq*($i/100))-1;
		print "$i\t$arr_read_count[$b_pi]\t$arr_minor_freq[$f_pi]\n";
}
print "max\t$arr_read_count[-1]\t$arr_minor_freq[-1]\n";
=cut
