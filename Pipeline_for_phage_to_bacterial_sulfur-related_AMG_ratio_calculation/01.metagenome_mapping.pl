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
my $adr = "/slowdata/data1/IMGVR/Transcriptomic_mapping/2nd_mapping";
#`cat /slowdata/data1/Lau/Lau_bins_NCBI_annotated/gbk_mdf_files/*.genome > All_genome_cat.genome`;
system ("bowtie2-build $adr/$metaG_id.a.gene.fasta $metaG_id.genome_scaffold");
open OUT, ">tmp_All_genome_cat_metagenome_mapping.txt";
foreach my $key (sort keys %map){
        if ($map{$key}[1] eq "DNA" and $map{$key}[3] eq $metaG_id){
			my $cmd = CMD($metaG_id, $key);
			print OUT $cmd;
		}
}
close OUT;

system ("cat tmp_All_genome_cat_metagenome_mapping.txt | parallel -j 50");
system ("rm tmp_All_genome_cat_metagenome_mapping.txt");
system ("rm $metaG_id.genome.*bt2");

}

sub CMD{
my $v1 = $_[0];
my $v2 = $_[1];
my $cmd = "bowtie2 -x ${v1}.genome_scaffold  -1 /slowdata/Reads/Lau/Plumes/original/${v2}_1.fastq  -2 /slowdata/Reads/Lau/Plumes/original/${v2}_2.fastq  -S ${v1}.genome_$map{$v2}[0]_mapped.sam -p 128;";
$cmd = $cmd."samtools view -bS ${v1}.genome_$map{$v2}[0]_mapped.sam > ${v1}.genome_$map{$v2}[0]_mapped.bam -@ 128;";
$cmd = $cmd."bamtools sort -in ${v1}.genome_$map{$v2}[0]_mapped.bam -out ${v1}.genome_$map{$v2}[0]_mapped.sorted.bam -mem 500;";
$cmd = $cmd."samtools index ${v1}.genome_$map{$v2}[0]_mapped.sorted.bam;";
$cmd = $cmd."samtools flagstat ${v1}.genome_$map{$v2}[0]_mapped.sorted.bam > ${v1}.genome_$map{$v2}[0]_mapped.sorted.stat;";
$cmd = $cmd."rm ${v1}.genome_$map{$v2}[0]_mapped.sam ${v1}.genome_$map{$v2}[0]_mapped.bam\n";
return $cmd;
}
