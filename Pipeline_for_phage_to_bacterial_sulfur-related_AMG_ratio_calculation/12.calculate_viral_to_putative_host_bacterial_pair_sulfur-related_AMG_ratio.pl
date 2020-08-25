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

my %Viral_Bacterial_pair = (); #pair =>  viral [0]  bacterial [1]
open IN, "paired_viral_and_bacterial_dsrA.txt";
while (<IN>){
	chomp;
	if (!/^#/){
		my @tmp = split (/\t/);
		$Viral_Bacterial_pair{$tmp[0]}[0] = $tmp[1];
		$Viral_Bacterial_pair{$tmp[0]}[1] = $tmp[2];
	}
}
close IN;

my %MetaG_ID_Map = (); # 3300001683 => Guaymas
open IN, "MetaG_ID_Map.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	$MetaG_ID_Map{$tmp[1]} = $tmp[0];
}
close IN;

foreach my $metaG_id (sort keys %MetaG){
	my %Depth = (); # gene => bam id => expression level
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
	
	my %Bam2pair2ratio = (); #bam file => pair =>  viral percentage, bacterial percentage, ratio result
	foreach my $bam (@Bam){
		my $viral_gene_coverage = 0;
		my $bacterial_gene_coverage = 0;
		foreach my $pair (sort keys %Viral_Bacterial_pair){
			if ($pair =~ /$MetaG_ID_Map{$metaG_id}/){
				my @Viral_genes = split (/\,/,$Viral_Bacterial_pair{$pair}[0]);
				my @Bacterial_genes = split (/\,/,$Viral_Bacterial_pair{$pair}[1]);
				my $pair_viral_gene_coverage = 0; my $pair_bacterial_gene_coverage = 0;
				my $pair_viral_gene_no = 0; my $pair_bacterial_gene_no = 0;
				foreach my $gene (sort keys %Gene_id){	
					if (exists $Viral_dsrA{$gene}){
						$viral_gene_coverage += $Depth{$gene}{$bam};
						if (@Viral_genes ~~ /$gene/){
							$pair_viral_gene_coverage += $Depth{$gene}{$bam}; $pair_viral_gene_no++;
						}
					}elsif (exists $Bacterial_dsrA{$gene}){
						$bacterial_gene_coverage += $Depth{$gene}{$bam};
						if (@Bacterial_genes ~~ /$gene/){
							$pair_bacterial_gene_coverage += $Depth{$gene}{$bam}; $pair_bacterial_gene_no++;
						}
					}
				
				}
			my $viral_percentage = $pair_viral_gene_coverage / (($viral_gene_coverage + $bacterial_gene_coverage) * $pair_viral_gene_no);
			my $bacterial_percentage = $pair_bacterial_gene_coverage / (($viral_gene_coverage + $bacterial_gene_coverage) * $pair_bacterial_gene_no);
			my $ratio = $viral_percentage / $bacterial_percentage;
			$Bam2pair2ratio{$bam}{$pair} = "$viral_percentage\,$bacterial_percentage\,$ratio";
				
			}# for each pair
		}
	}
	
	open OUT, ">$metaG_id.viral_to_bacterial.pair_info.dsrA_gene_ratio_result.txt";
	my $row=join("\t", @Bam);
	print OUT "Head\t$row\n";
	foreach my $tmp1 (sort keys %Viral_Bacterial_pair)	
	{
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (@Bam)
        {
                if (exists $Bam2pair2ratio{$tmp2}{$tmp1})
                {
                        push @tmp, $Bam2pair2ratio{$tmp2}{$tmp1};
                }
                else
                {
                        push @tmp,"0"
                }
        }
        print OUT join("\t",@tmp)."\n";
	}
	close OUT;
}


