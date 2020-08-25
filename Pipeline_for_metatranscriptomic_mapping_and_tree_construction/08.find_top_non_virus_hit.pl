#!/usr/bin/perl

use strict;
use warnings;

my %Unique = (); # acc => $anno [0], $tax [1]
open IN, "dsrA_total.aa.fasta_blast_result_col2_unique_grep_nr_remove_SOH";
while(<IN>){
	chomp;
	my $line = $_;
	my $acc = my $anno = my $tax = "";
	if($line =~ /MULTISPECIES\:/){
		($acc,$anno,$tax) = $line =~ /^>(.+?)\sMULTISPECIES\:\s(.+?)\s\[(.+)\]$/;  
	}elsif($line =~ /RecName/){
		$tax = "Unknown";
		if ($line !~ /\;/){
			($acc,$anno) = $line =~ /^>(.+?)\sRecName\:\sFull\=(.+)$/; 
		}else{
			($acc,$anno) = $line =~ /^>(.+?)\sRecName\:\sFull\=(.+?)\;/; 
		}
	}elsif($line =~ / Chain /){
		$tax = "Unknown";
		($acc) = $line =~ /^>(.+?)\sChain/;  
		($anno) = $line =~ /\,\s(.+)$/;  
	}elsif($line =~ /\|\|/){
		$tax = "Unknown";
		($acc,$anno) = $line =~ /^>(.+?)\s(.+)$/;  
	}elsif($line =~ /\.\d\s\[/){
		($acc,$tax) = $line =~ /^>(.+?)\s\[(.+)\]$/; $anno = "hypothetical protein";
	}else{
		($acc,$anno,$tax) = $line =~ /^>(.+?)\s(.+?)\s\[(.+)\]$/;
	}
	$Unique{$acc}[0] = $anno;$Unique{$acc}[1]  = $tax;
}
close IN;

`mkdir tmp`;
my $Individual_hits = (); #$gene => $line
my $gene = ""; 
open IN, "dsrA_total.aa.fasta_blast_result.txt";
while (<IN>){
	chomp;
	my $line = $_;
	my @tmp = split (/\t/);
	if (!$gene){
		$gene = $tmp[0];
		$Individual_hits .= "$line";
	}elsif($gene eq $tmp[0]){
		$Individual_hits .= "\n$line";
	}else{
		open OUT, ">tmp/dsrA_total.aa.fasta_blast_result_Subfile_$gene.txt";
		print OUT "$Individual_hits\n";
		close OUT;
		$gene = ""; $Individual_hits = "";
	}
}
close IN;

my %All_top_hits = (); #gene => [0] $acc, [1] $iden, [2] $e, [3] $bitscore  non-virus hits
open INN, "ls tmp/*.txt|";
while (<INN>){
	chomp;
	my $file = $_;
	my ($gene) = $_ =~ /^.+?Subfile\_(.+?)\.txt/;
	my $acc = my $iden = my $e = my $bitscore = ""; #the non-virus top hit 
	my $last_line = "";
	open IN, "$file";
	while (<IN>){
		chomp;
		my @tmp = split (/\t/);
		if ($Unique{$tmp[1]}[1] !~ /virus/){
			$acc = $tmp[1];
			$iden = $tmp[2];
			$e = $tmp[-2];
			$bitscore = $tmp[-1];
			last;
		}
		$last_line = $_ if eof;
	}
	close IN;
	if (!$acc){
		my @tmp = split (/\t/, $last_line);
		$acc = $tmp[1];
		$iden = $tmp[2];
		$e = $tmp[-2];
		$bitscore = $tmp[-1];
	}
	$All_top_hits{$gene}[0] = $acc;
	$All_top_hits{$gene}[1] = $iden;
	$All_top_hits{$gene}[2] = $e;
	$All_top_hits{$gene}[3] = $bitscore;
}
close INN;


open OUT, ">dsrA_total.aa.fasta_blast_result_parsed.txt_blastp_new_nr_non-virus_tophit.txt";
print OUT "query\taccession\tproduct\ttaxonomy\tidentity\tevalue\tbitscore\n";
foreach my $gene (sort keys %All_top_hits){
	my $acc = $All_top_hits{$gene}[0];
	my $iden = $All_top_hits{$gene}[1];
	my $e = $All_top_hits{$gene}[2];
	my $bitscore = $All_top_hits{$gene}[3];
	my $anno = my $tax = "";
	$anno = $Unique{$acc}[0] or $anno = "NA";
	$tax = $Unique{$acc}[1] or $tax = "NA";
	print OUT "$gene\t$acc\t$anno\t$tax\t$iden\t$e\t$bitscore\n";
}
close OUT;
