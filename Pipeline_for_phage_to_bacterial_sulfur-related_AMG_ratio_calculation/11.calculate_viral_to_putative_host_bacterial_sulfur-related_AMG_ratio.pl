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

my %Viral_dsrA_SUP05_like_Clade1 = (); #gene id => 1
my %Bacterial_dsrA_SUP05_like_Clade1 = ();#gene id => 1
open IN, "dsrA_gene_phylogeny.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	if ($tmp[0] =~ /__bacterial_dsrA/ and $tmp[1] eq "SUP05_like-Clade1"){
		my ($gene_id) = $tmp[0] =~ /^(.+?)__bacterial_dsrA/;
		$Bacterial_dsrA_SUP05_like_Clade1{$gene_id} = 1; 
	}elsif($tmp[0] =~ /__viral_dsrA/ and $tmp[1] eq "SUP05_like-Clade1"){
		my ($gene_id) = $tmp[0] =~ /^(.+?)__viral_dsrA/;
		$Viral_dsrA_SUP05_like_Clade1{$gene_id} = 1;
	}
}
close IN;

my %Viral_dsrA_SUP05_like_Clade2 = (); #gene id => 1
my %Bacterial_dsrA_SUP05_like_Clade2 = ();#gene id => 1
open IN, "dsrA_gene_phylogeny.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	if ($tmp[0] =~ /__bacterial_dsrA/ and $tmp[1] eq "SUP05_like-Clade2"){
		my ($gene_id) = $tmp[0] =~ /^(.+?)__bacterial_dsrA/;
		$Bacterial_dsrA_SUP05_like_Clade2{$gene_id} = 1;
	}elsif($tmp[0] =~ /__viral_dsrA/ and $tmp[1] eq "SUP05_like-Clade2"){
		my ($gene_id) = $tmp[0] =~ /^(.+?)__viral_dsrA/;
		$Viral_dsrA_SUP05_like_Clade2{$gene_id} = 1;
	}
}
close IN;


foreach my $metaG_id (sort keys %MetaG){
	my %Depth = (); # gene => bam file => expression level
	my @Head = (); my @Bam = (); my %Gene_id = ();
	open IN, "$metaG_id.genome.depth.txt";
	while (<IN>){
		chomp;
		if (/^contigName/){
			@Head = split (/\t/);
			foreach my $key (@Head){
				if ($key =~ /\.sorted\.bam$/){
					push @Bam,$key;
				}
			}
		}else{
			my @tmp = split (/\t/);
			for(my $i=1; $i<=$#tmp; $i++){
				$Depth{$tmp[0]}{$Head[$i]} = $tmp[$i]; $Gene_id{$tmp[0]} = 1;
			}
		}
	}
	close IN;
	
	my %Bam2ratio = (); #bam file => ratio result
	foreach my $bam (@Bam){
		my $viral_gene_coverage = 0;
		my $bacterial_gene_coverage = 0;
		my $viral_gene_coverage_Clade1 = 0;
		my $bacterial_gene_coverage_Clade1 = 0;
		my $viral_gene_coverage_Clade2 = 0;
		my $bacterial_gene_coverage_Clade2 = 0;
		foreach my $gene (sort keys %Gene_id){					
			if (exists $Viral_dsrA{$gene}){
				$viral_gene_coverage += $Depth{$gene}{$bam};
				if (exists $Viral_dsrA_SUP05_like_Clade1{$gene}){
					$viral_gene_coverage_Clade1 += $Depth{$gene}{$bam};
				}elsif (exists $Viral_dsrA_SUP05_like_Clade2{$gene}){
					$viral_gene_coverage_Clade2 += $Depth{$gene}{$bam};
				}
			}		
			if (exists $Bacterial_dsrA{$gene}){
				$bacterial_gene_coverage += $Depth{$gene}{$bam};
				if (exists $Bacterial_dsrA_SUP05_like_Clade1{$gene}){
					$bacterial_gene_coverage_Clade1 += $Depth{$gene}{$bam};
				}elsif ($Bacterial_dsrA_SUP05_like_Clade2{$gene}){
					$bacterial_gene_coverage_Clade2 += $Depth{$gene}{$bam};
				}
			}
		}
		my $ratio1 = $viral_gene_coverage_Clade1/($viral_gene_coverage + $bacterial_gene_coverage);
		my $ratio2 = $viral_gene_coverage_Clade2/($viral_gene_coverage + $bacterial_gene_coverage);
		my $ratio3 = $bacterial_gene_coverage_Clade1/($viral_gene_coverage + $bacterial_gene_coverage);
		my $ratio4 = $bacterial_gene_coverage_Clade2/($viral_gene_coverage + $bacterial_gene_coverage);
		$Bam2ratio{$bam} = "$ratio1\t$ratio2\t$ratio3\t$ratio4";
	}
	
	open OUT, ">$metaG_id.viral_to_bacterial.clade_info.dsrA_gene_ratio_result.txt";
	print OUT "\tviral_Clade1_percentage\tviral_Clade2_percentage\tbacterial_Clade1_percentage\tbacterial_Clade2_percentage\n";
	foreach my $key (sort keys %Bam2ratio){
		print OUT "$key\t$Bam2ratio{$key}\n";
	}
	close OUT;	
}


