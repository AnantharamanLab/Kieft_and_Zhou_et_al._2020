for i in `ls All_gene_cat.gene*mapped.sorted.bam`
do
/slowdata/archive/bbmap/pileup.sh in=$i out=$i.pileup.out
done
