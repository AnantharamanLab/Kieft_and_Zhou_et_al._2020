#!/usr/bin/perl -w

open INN, "ls *pileup.out | ";
while (<INN>){
	chomp; my $pileup_out = $_;
	open IN, "$pileup_out";
	my %hash = ();
	while (<IN>)
	{
		chomp;
		unless (/#/)
		{
			my @tmp = split(/\t/, $_);
			$hash{$tmp[0]}[0] = $tmp[6] + $tmp[7];
			$hash{$tmp[0]}[1] = $tmp[2];	
		}
	}
	close IN;

	open OUT, ">$pileup_out.parsed.csv";
	foreach my $key (sort keys %hash)
	{
		print OUT "$key\,$hash{$key}[0]\,$hash{$key}[1]\n";
	}
	close OUT;
} # to open each pileup files
close INN;
