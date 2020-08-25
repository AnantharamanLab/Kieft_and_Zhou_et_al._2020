# Kieft_and_Zhou_et_al._2020

This GitHub repository contains supplementary files associated with the manuscript _Ecology of inorganic sulfur auxiliary metabolism in widespread bacteriophages_

The bioRxiv preprint of this manuscript can be found [here](https://www.biorxiv.org/content/10.1101/2020.08.24.253096v1). 

8/25/2020  
Zhichao Zhou  
Kristopher Kieft  
Karthik Anantharaman  
University of Wisconsin-Madison  


## Table of Contents:
1. [Data files](#data)
2. [Taxonomy](#taxonomy)
3. [Phage to bacterial sulfur-related AMG ratio calculation](#map)
4. [Metatranscriptomic mapping](#ratio)
5. [Contact](#contact)


## Explanations


### Data Files <a name="data"></a>

* `Kieft_and_Zhou_et_al_2020.genomes.fasta`: all 191 vMAG genome sequences used in this study  
* `Kieft_and_Zhou_et_al_2020.proteins.faa`: all predicted protein sequences for the vMAGs used in this study  


### Taxonomy <a name="taxonomy"></a>
This folder contains the custom Python program and associated database to estimate the taxonomy of prokaryotic viruses. The program (`PTT.py`, Phage Taxonomy Tool) takes either proteins or genomes as an input file and will run Diamond Blastp against a database of reference viruses from NCBI GenBank. Diamond Blastp hits will be filtered and assignment of a taxonomic lineage will be determined according to counting Diamond Blastp hits. Assignments will be made on the levels of Order, Family and Sub-family, though database metadata is provided to the Genus level if manual curation of Genus is desired. Each taxonomic level will be assigned individually, starting with Order. If an Order is assigned it will move to Family before assigning Sub-family. For example, if Order _Caudovirales_ and Family _unknown_ is assigned then there will be no assignment of Sub-family since the higher level (i.e., Family) was unknown.  

_Dependencies_: Diamond, Prodigal (if input is nucleotides), Blast, BioPython. 

_Before starting_: Please unzip the Diamond database file using the command `gunzip PTT_database.dmnd.gz`. In the event you are using an uncompatible version of Diamond (v2 used here) you may need to generate the database on your own. Please unzip the protein database using the command `gunzip PTT_database.faa.gz` followed by the command `diamond makedb --in PTT_database.faa --db PTT_database.dmnd`. In addition to the Diamond database there is also a taxonomy metadata file required for the program (`PTT_virus_taxonomy.tsv`). Please place the database files (`PTT_database.dmnd` and `PTT_virus_taxonomy.tsv`) in the same location as the program script (`PTT.py`). If the databases are not in the same location you may choose to use the optional flags `-d` and `-v` to designate the database location(s).  

_Usage_: This programs takes a minimum of one argument (the input file). For quick usage, run `python3 PTT.py -i input_file.faa` with `input_file.faa` being a protein fasta file of your choice in Prodigal format. `-t` to increase the number of threads used. `-f` to indicate "prot" (protein) or "nucl" (nucleotide) input. `-d` to provide the location for the Diamond database (if not in same folder as `PTT.py`). `-v` to provide the location for the metadata taxonomy file (if not in same folder as `PTT.py`). See `PTT.py -h` for the help menu. Note, Diamond is very fast and in most cases a single thread is sufficient runtime. For example, runtime of all proteins from this manuscript with `-t 1` is approximately 25 seconds.  

_Outputs_: An automatic output folder will be generated named `PTT_input_date` where input comes from the `-i` flag and date is the date of program finish. The output folder will contain one important file, `input.PTT.virus-taxonomy.tsv`. This will contain tab separated Order, Family and Sub-family assignments for each input genome. The column "pred_score" contains 3 numbers, each corresponding to the Order/Family/Sub-family predictions and have a maximum value of 6, but may be ignored. The file `input.PTT.protein-taxonomy.tsv` contains the same information for each protein individually, with the addition of a Genus prediction. The log file contains general run information. The folder `analysis_files` contains general information used for analysis including raw/parsed Diamond ouputs. The `tophits.tsv` and `toppreps.tsv` files may be ignored.  

_Output Terminology_: Several assignment terms may be seen in addition to standard nomenclature (standard nomenclature examples: _Caudovirales_, _Myoviridae_). These terms include _unassigned_, _ambiguous_ and _unknown_. The term _unassigned_ can be found in both the output files and `PTT_virus_taxonomy.tsv` metadata. This term refers to taxonomic levels that were not assigned according to the NCBI GenBank database. Therefore, the program may have identified a significant hit to a particular protein or virus, but the database contains "unassigned" information. The term _ambiguous_ referes to a situation in which a protein or virus had significant hits to multiple proteins within the database, but the program could not distinguish between taxonomic assignments. For example, this occurs if the program cannot distinguish an input virus as belonging to _Myoviridae_ or _Siphoviridae_. In these cases manual verification can be used. The final term is _unknown_ which indicates situtations in which there were not a significant number of hits to the protein database and a scaffold was unable to be assigned to a taxonomic level nor be assigned _ambiguous_. 

_Cautions_: The input file headers (sequence names) cannot have spaces in the names or else the output will cut the names at those spaces. Input proteins must be in Prodigal format (e.g., >sequence_name_#). The program can only assign taxonomy that is contained within the protein database and taxonomy metadata files. Therefore, it is unsuitable for all eukaryotic taxonomic assignments. It is also unsuitable for novel phage Families, such as _Chaseviridae_ within the Order _Caudovirales_. 



### Phage to bacterial sulfur-related AMG ratio calculation <a name="map"></a>

_Dependencies_: Perl v5+, Bowtie 2 v2.3.4.1,  jgi_summarize_bam_contig_depths (implemented in MetaWrap), Diamond v0.9.28.129, MAFFT v7.271, IQ-TREE v1.6.9

_Explanation to each step in the pipeline:_    
01.metagenome_mapping.pl    
Concatenate all the genes from metagenome as the mapping reference; Use Bowtie 2 to map filterd and QC-processed reads; Finally get sorted bam files as the result. The metadata file "Metagenome_map.txt" was used to allow processing mutiple metagenomes mapping by this script.    

02.calculate_metagenome_to_MAG_depth.pl    
Use "jgi_summarize_bam_contig_depths" to calculate gene coverage.    

03.grep_dsrA_list.sh    
Grep _dsrA_ genes from metagenome. Each IMG metagenome has its annotation by the DOE IMG database. We used the annotation to pre-select _dsrA_ genes from each metagenome. It needs further manual curation.  

04.calculate_viral_to_bacterial_sulfur-related_AMG_ratio.pl    
Read gene coverage result from Step #2. Calculate the viral to bacterial total sulfur-related AMG gene coverage ratio for each metagenome.    

05.grep_all_dsrA_gene.pl    
Grep all the dsrA gene encoding proteins from the metagenome.    

06.run_blastp.sh    
Run Diamond Blastp for all the DsrA sequences.    

07.parse_blast_result.sh    
Parse the Diamond Blastp result.    

08.find_top_non_virus_hit.pl    
Parse the Diamond Blastp result to screen for top 10 non-virus hits, which will be downloaded and used as reference sequences to build phylogenetic tree.    

09.make_dsrA_tree.pl    
Use viral and bacterial DsrA sequences, and reference DsrA sequences from Step #8, to build phylogenetic tree. MAFFT was used to align the sequences, and IQ-TREE was used to build the tree.    

10.Replace_tip_names.pl    
Replace tip names of resulted tree to the ones that are meaningful and formal.       

11.calculate_viral_to_bacterial_sulfur-related_AMG_ratio_according_to_taxonomy.pl
Firstly, divide the viral and bacterial sequences to each category according to their taxonomy. Then, calculate their AMG gene coverage ratios.    

12.calculate_viral_to_putative_host_bacterial_pair_sulfur-related_AMG_ratio.pl        
Firstly, get the viral to putative host bacterial pair information. Then, calculate viral and bacterial sequence gene coverage percentage values within each pair (the gene coverage values of viral and bacterial sequence were all normalized by gene numbers). Fianlly, calcualte the viral to putative host bacterial gene coverage ratios.


### Metatranscriptomic mapping <a name="ratio"></a>

_Dependencies_: Perl v5+, Bowtie 2 v2.3.4.1,  pileup (implemented in BBmap), Diamond v0.9.28.129, MAFFT v7.271, IQ-TREE v1.6.9    

01.Transcriptom_mapping.pl    
Concatenate all the genes from metagenome as the mapping reference; Use Bowtie 2 to map filterd, QC-processed and rRNA removed metatranscriptomic reads; Finally get sorted bam files as the result. The metadata file "Transcriptome_map.txt" was used to allow processing mutiple metatranscriptomes mapping by this script.        

02.pileup_to_calculate_gene_reads_abundance_for_MetaT.sh
Use "pileup.sh" to calculate gene abundance.   

03.parse_pileup_info.pl    
Parse pileup result.    

04.calculate_MetaT_RPKM.pl    
Calculate the gene expression result by normalizing gene abundance by read numbers (1 million reads) and gene length (1k).    

05.parse_MetaT_RPKM_result.pl    
Parse to get MetaT RPKM results for targeted genes.    

06.run_blastp.sh    
Run Diamond Blastp for all the DsrA sequences.    

07.parse_blast_result.sh    
Parse the Diamond Blastp result.    

08.find_top_non_virus_hit.pl     
Parse the Diamond Blastp result to screen for top 10 non-virus hits, which will be downloaded and used as reference sequences to build phylogenetic tree.     

09.make_dsrA_tree.pl    
Use viral and bacterial DsrA sequences, and reference DsrA sequences from Step #8, to build phylogenetic tree. MAFFT was used to align the sequences, and IQ-TREE was used to build the tree.     

10.Replace_tip_names.pl    
Replace tip names of resulted tree to the ones that are meaningful and formal.    

### Contact <a name="contact"></a>

Kristopher Kieft, kieft@wisc.edu  
Zhichao Zhou, zzhou388@wisc.edu 
