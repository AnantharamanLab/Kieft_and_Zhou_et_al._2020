#!/usr/bin/perl

use strict;
use warnings;

my %Target_gene_list = ();
open IN, "gene.list";
while (<IN>){
	chomp;
	$Target_gene_list{$_} = 1;
}
close IN;

my $head = ""; my %MetaT_RPKM = ();
open IN, "MetaT.RPKM.txt";
while (<IN>){
	chomp;
	if (/^Head/){
		$head = $_;
	}else{
		my @tmp = split (/\t/);
		$MetaT_RPKM{$tmp[0]} = $_;
	}
}
close IN;

open OUT, ">MetaT_RPKM-target_gene.txt";
print OUT $head."\n";
foreach my $key (sort keys %MetaT_RPKM){
	if (exists $Target_gene_list{$key}){
		print OUT $MetaT_RPKM{$key}."\n";
	}
}
close OUT;
