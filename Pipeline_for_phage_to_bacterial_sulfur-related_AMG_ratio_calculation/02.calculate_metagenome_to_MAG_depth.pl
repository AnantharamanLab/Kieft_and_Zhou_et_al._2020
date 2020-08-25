#!/usr/bin/perl

use strict;
use warnings;
my %MetaG = (); 
my %map = (); # DNA and cDNA reads maps
open IN, "Metagenome_map.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	$map{$tmp[3]}[0] = $tmp[0];   # SRR [0] => name
	$map{$tmp[3]}[1] = $tmp[1];   # SRR [1] => cDNA or DNA
	$map{$tmp[3]}[2] = $tmp[2];   # SRR [2] => PAIRED or SINGLE
	$map{$tmp[3]}[3] = $tmp[-1];   # SRR [3] => corresponding assembly
	$MetaG{$tmp[-1]} = 1;
}
close IN;

foreach my $metaG_id (sort keys %MetaG){
	`jgi_summarize_bam_contig_depths --outputDepth $metaG_id.genome.depth.txt  --pairedContigs $metaG_id.genome.paired.txt  $metaG_id.*.sorted.bam`;
}
