#!/usr/bin/perl

use strict;
use warnings; 

my %map = (); # DNA and cDNA reads maps
open IN, "Transcriptome_map.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	$map{$tmp[3]}[0] = $tmp[0];   # SRR [0] => name
	$map{$tmp[3]}[1] = $tmp[1];   # SRR [1] => cDNA or DNA
	$map{$tmp[3]}[2] = $tmp[2];   # SRR [2] => PAIRED or SINGLE
	$map{$tmp[3]}[3] = $tmp[4];   # SRR [3] => reads numbers
}
close IN;

my $gn = "All_gene_cat"; 
#`cp /slowdata/data2/IMG_Phages/Dissimilatory/IMGVR_sox_dsr/Genomes_proteins/IMGVR_sulfur_AMGs_phages_genes-meta.fna All_gene_cat.gene`;
system ("bowtie2-build $gn.gene $gn.gene_scaffold");
open OUT, ">tmp_All_gene_cat_transcriptom_mapping.txt";
foreach my $key (sort keys %map){
        if ($map{$key}[1] eq "cDNA" and $map{$key}[2] eq "SINGLE"){
			my $cmd = CMD1($gn, $key);
			print OUT $cmd;
		}elsif ($map{$key}[1] eq "cDNA" and $map{$key}[2] eq "PAIRED"){
			my $cmd = CMD2($gn, $key);
			print OUT $cmd;
		}
}
close OUT;

system ("cat tmp_All_gene_cat_transcriptom_mapping.txt | parallel -j 50");
system ("rm tmp_All_gene_cat_transcriptom_mapping.txt");
#system ("rm $gn.gene_scaffold.*bt2");

#add parameters to let bowtie2 only find properly mapped reads (excluding single read mapping, unproperly mapping)
#and use the loose criterion to find mapping reads (which will report more potential reads)
#did not use the local mode, due to that I am afraid to find reads that is mapping to the boundary of genes and surpass the end of genes. (This point I am not sure about it)
#More information : Chinese explantation on parameters of bowtie2 : https://www.plob.org/article/4540.html
#English manual of details on the use of bowtie2 : http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#end-to-end-alignment-versus-local-alignment

#in the 1st test, I use "--no-mixed --no-discordant  --no-dovetail --no-contain --no-overlap --very-sensitive" to find mapped reads
#in the 2nd test, I use "--very-sensitive" to find mapped reads
sub CMD1{
my $v1 = $_[0];
my $v2 = $_[1];
my $cmd = "bowtie2 -x ${v1}.gene_scaffold  -U /slowdata/Reads/Lau/MetaT_reads/${v2}.trim_non_rRNA.fastq  -S ${v1}.gene_$map{$v2}[0]_mapped.sam -p 128  --very-sensitive;";
$cmd = $cmd."samtools view -bS ${v1}.gene_$map{$v2}[0]_mapped.sam > ${v1}.gene_$map{$v2}[0]_mapped.bam -@ 128;";
$cmd = $cmd."bamtools sort -in ${v1}.gene_$map{$v2}[0]_mapped.bam -out ${v1}.gene_$map{$v2}[0]_mapped.sorted.bam -mem 500;";
$cmd = $cmd."samtools index ${v1}.gene_$map{$v2}[0]_mapped.sorted.bam;";
$cmd = $cmd."samtools flagstat ${v1}.gene_$map{$v2}[0]_mapped.sorted.bam > ${v1}.gene_$map{$v2}[0]_mapped.sorted.stat;";
$cmd = $cmd."rm ${v1}.gene_$map{$v2}[0]_mapped.sam ${v1}.gene_$map{$v2}[0]_mapped.bam\n";
return $cmd;
}

sub CMD2{
my $v1 = $_[0];
my $v2 = $_[1];
my $cmd = "bowtie2 -x ${v1}.gene_scaffold  -1 /slowdata/Reads/Lau/MetaT_reads/${v2}.trim_non_rRNA.R1.fastq  -2 /slowdata/Reads/Lau/MetaT_reads/${v2}.trim_non_rRNA.R2.fastq -S ${v1}.gene_$map{$v2}[0]_mapped.sam -p 128  --very-sensitive;";
$cmd = $cmd."samtools view -bS ${v1}.gene_$map{$v2}[0]_mapped.sam > ${v1}.gene_$map{$v2}[0]_mapped.bam -@ 128;";
$cmd = $cmd."bamtools sort -in ${v1}.gene_$map{$v2}[0]_mapped.bam -out ${v1}.gene_$map{$v2}[0]_mapped.sorted.bam -mem 500;";
$cmd = $cmd."samtools index ${v1}.gene_$map{$v2}[0]_mapped.sorted.bam;";
$cmd = $cmd."samtools flagstat ${v1}.gene_$map{$v2}[0]_mapped.sorted.bam > ${v1}.gene_$map{$v2}[0]_mapped.sorted.stat;";
$cmd = $cmd."rm ${v1}.gene_$map{$v2}[0]_mapped.sam ${v1}.gene_$map{$v2}[0]_mapped.bam\n";
return $cmd;
}

