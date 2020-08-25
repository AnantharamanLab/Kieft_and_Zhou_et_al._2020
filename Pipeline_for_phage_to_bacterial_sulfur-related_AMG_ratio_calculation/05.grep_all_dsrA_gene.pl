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

my %Seq = (); my $head = ();
foreach my $metaG_id (sort keys %MetaG){
	open IN, "$metaG_id/$metaG_id.a.faa";
	while (<IN>){
		chomp;
		if (/>/){
			$head = $_;
			$Seq{$head} = "";
		}else{
			$Seq{$head} .= $_;
		}
	}
	close IN;
}

my %Viral_dsrA = (); #gene id => metagenome ID
open IN, "Lau_viral_sulfur_related_AMG_list.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	$Viral_dsrA{$tmp[2]} = $tmp[0]; #print "$tmp[2]\n";
}
close IN;

my %Bacterial_dsrA = ();#gene id => metagenome ID
open IN, "All_dsrA_gene_list.txt";
while (<IN>){
	chomp;
	if (!exists $Viral_dsrA{$_}){
		$Bacterial_dsrA{$_} = 1;
	}
}
close IN;

open OUT, ">Lau_dsrA_protein.faa";
foreach my $key (sort keys %Seq){
	my ($gene_id) = $key =~ />(.+)$/;
	if (exists $Viral_dsrA{$gene_id}){
		print OUT "$key"."__viral_dsrA\n";
		my $fasta_seq = $Seq{$key};
		$fasta_seq =~ s/\*$//g;
		print OUT $fasta_seq."\n";
	}elsif (exists $Bacterial_dsrA{$gene_id}){
		print OUT "$key"."__bacterial_dsrA\n";
		my $fasta_seq = $Seq{$key};
		$fasta_seq =~ s/\*$//g;
		print OUT $fasta_seq."\n";
	}
}
close OUT;



