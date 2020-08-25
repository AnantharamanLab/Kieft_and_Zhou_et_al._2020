#!/usr/bin/perl

use strict;
use warnings;

#1."per million" scaling factor is calculated as total number of mapped reads divided by 1,000,000
#2. Divide the read counts by the "per million" scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM)
#3. Divide the RPM values by the length of the gene, in kilobases. This gives you RPKM.

my %map = (); # DNA and cDNA reads maps
open IN, "Transcriptom_map.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	$map{$tmp[3]}[0] = $tmp[0];   # SRR [0] => name
	$map{$tmp[3]}[1] = $tmp[1];   # SRR [1] => cDNA or DNA
	$map{$tmp[3]}[2] = $tmp[2];   # SRR [2] => PAIRED or SINGLE
	$map{$tmp[3]}[3] = $tmp[4];   # SRR [3] => reads numbers
}
close IN;

my %h = (); my %MetaT = ();
open IN, "ls All_gene_cat.gene_*_mapped.sorted.bam.pileup.out.parsed.csv | ";
while (<IN>){
	chomp;
	my $file = $_;
	my ($name) = $_ =~ /^All\_gene\_cat\.gene\_(.+?)\_mapped\.sorted\.bam\.pileup\.out\.parsed\.csv$/;
	$MetaT{$name} =1;
	my $reads_num = 0;
	foreach my $key (sort keys %map){
		if ($map{$key}[0] eq $name){
			$reads_num = $map{$key}[3];
		}
	}
	
	open INN, $file;
	while (<INN>){
		chomp;
		my @tmp = split (/\,/,$_);
		my $gene_read_count = $tmp[1]; my $gene_read_len = $tmp[2]; 
		$h{$tmp[0]}{$name} = $gene_read_count * (1000000 / $reads_num) / ($gene_read_len /1000);
	}
	close INN;

}
close IN;
	
#print table
open OUT, ">MetaT.RPKM.txt";
my $row=join("\t", sort keys %MetaT);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %h)
{
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %MetaT)
        {
                if (exists $h{$tmp1}{$tmp2})
                {
                        push @tmp, $h{$tmp1}{$tmp2};
                }
                else
                {
                        push @tmp,"0"
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;
