i="dsrA_total.aa.fasta"

cut -f2 ${i}_blast_result.txt | sort -u > ${i}_blast_result_col2_unique
grep -Ff ${i}_blast_result_col2_unique /slowdata/databases/NCBI_nr_diamond/nr > ${i}_blast_result_col2_unique_grep_nr
cp ${i}_blast_result_col2_unique_grep_nr ${i}_blast_result_col2_unique_grep_nr_remove_SOH
sed -i "s/\x01.*$//" ${i}_blast_result_col2_unique_grep_nr_remove_SOH
#Rscript --vanilla thor_3_bins_blastp_new_nr.R ${i}_blast_result.txt ${i}_blast_result_col2_unique_grep_nr_remove_SOH ${i}_blast_result_parsed.txt

