#!/usr/bin/perl

use strict;
use warnings;

my %hash = (); # UWM id => strain name
open IN, "tip_map.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/,$_);
	$hash{$tmp[0]} = $tmp[1];
}
close IN;

open IN, "Lau_dsrA_protein_and_ref.mafft.faa.treefile";
open OUT, ">Lau_dsrA_protein_and_ref.mafft.faa.treefile.tip_name_mdf.nwk";
while (<IN>){
	chomp;
	foreach my $key (sort keys %hash){
		if ($_ =~ /$key/){
			$_ =~ s/$key\:/$hash{$key}\:/g;
		}
	}
	print OUT "$_\n";
}
close IN;
close OUT;

