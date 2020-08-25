# Kieft_and_Zhou_et_al._2020

This GitHub repository contains supplementary files associated with the manuscript _Ecology of inorganic sulfur auxiliary metabolism in widespread bacteriophages_

The bioRxiv preprint of this manuscript can be found [here](https://www.biorxiv.org/content/10.1101/2020.08.24.253096v1). 

8/25/2020  
Zhichao Zhou  
Kristopher Kieft  
Karthik Anantharaman  
University of Wisconsin-Madison  


## Table of Contents:
1. [Taxonomy](#taxonomy)
2. [Mapping and Tree Construction](#map)
3. [Phage to Bacteria Ratio Calculation](#ratio)
4. [Contact](#contact)


## Explanations

### Taxonomy <a name="taxonomy"></a>
This folder contains the custom Python program and associated database to estimate the taxonomy of prokaryotic viruses. The program (`PTT.py`, Phage Taxonomy Tool) takes either proteins or genomes as an input file and will run Diamond Blastp against a database of reference viruses from NCBI GenBank. Diamond Blastp hits will be filtered and assignment of a taxonomic lineage will be determined according to counting Diamond Blastp hits. Assignments will be made on the levels of Order, Family and Sub-family, though database metadata is provided to the Genus level if manual curation of Genus is desired. Each taxonomic level will be assigned individually, starting with Order. If an Order is assigned it will move to Family before assigning Sub-family. For example, if Order _Caudovirales_ and Family _unknown_ is assigned then there will be no assignment of Sub-family since the higher level (i.e., Family) was unknown.  

_Dependencies_: Diamond, Prodigal (if input is nucleotides), Blast, BioPython. 

_Before starting_: Please unzip the Diamond database file using the command `gunzip PTT_database.dmnd.gz`. In the event you are using an uncompatible version of Diamond (v2 used here) you may need to generate the database on your own. Please unzip the protein database using the command `gunzip PTT_database.faa.gz` followed by the command `diamond makedb --in PTT_database.faa --db PTT_database.dmnd`. In addition to the Diamond database there is also a taxonomy metadata file required for the program (`PTT_virus_taxonomy.tsv`). Please place the database files (`PTT_database.dmnd` and `PTT_virus_taxonomy.tsv`) in the same location as the program script (`PTT.py`). If the databases are not in the same location you may choose to use the optional flags `-d` and `-v` to designate the database location(s).  

_Usage_: This programs takes a minimum of one argument (the input file). For quick usage, run `python3 PTT.py -i input_file.faa` with `input_file.faa` being a protein fasta file of your choice in Prodigal format. `-t` to increase the number of threads used. `-f` to indicate "prot" (protein) or "nucl" (nucleotide) input. `-d` to provide the location for the Diamond database (if not in same folder as `PTT.py`). `-v` to provide the location for the metadata taxonomy file (if not in same folder as `PTT.py`). See `PTT.py -h` for the help menu. Note, Diamond is very fast and in most cases a single thread is sufficient runtime. For example, runtime of all proteins from this manuscript with `-t 1` is approximately 25 seconds.  

_Outputs_: An automatic output folder will be generated named `PTT_input_date` where input comes from the `-i` flag and date is the date of program finish. The output folder will contain one important file, `input.PTT.virus-taxonomy.tsv`. This will contain tab separated Order, Family and Sub-family assignments for each input genome. The column "pred_score" contains 3 numbers, each corresponding to the Order/Family/Sub-family predictions and have a maximum value of 6, but may be ignored. The file `input.PTT.protein-taxonomy.tsv` contains the same information for each protein individually, with the addition of a Genus prediction. The log file contains general run information. The folder `analysis_files` contains general information used for analysis including raw/parsed Diamond ouputs. The `tophits.tsv` and `toppreps.tsv` files may be ignored.  

_Output Terminology_: Several assignment terms may be seen in addition to standard nomenclature (standard nomenclature examples: _Caudovirales_, _Myoviridae_). These terms include _unassigned_, _ambiguous_ and _unknown_. The term _unassigned_ can be found in both the output files and `PTT_virus_taxonomy.tsv` metadata. This term refers to taxonomic levels that were not assigned according to the NCBI GenBank database. Therefore, the program may have identified a significant hit to a particular protein or virus, but the database contains "unassigned" information. The term _ambiguous_ referes to a situation in which a protein or virus had significant hits to multiple proteins within the database, but the program could not distinguish between taxonomic assignments. For example, this occurs if the program cannot distinguish an input virus as belonging to _Myoviridae_ or _Siphoviridae_. In these cases manual verification can be used. The final term is _unknown_ which indicates situtations in which there were not a significant number of hits to the protein database and a scaffold was unable to be assigned to a taxonomic level nor be assigned _ambiguous_. 

_Cautions_: The input file headers (sequence names) cannot have spaces in the names or else the output will cut the names at those spaces. Input proteins must be in Prodigal format (e.g., >sequence_name_#). The program can only assign taxonomy that is contained within the protein database and taxonomy metadata files. Therefore, it is unsuitable for all eukaryotic taxonomic assignments. It is also unsuitable for novel phage Families, such as _Chaseviridae_ within the Order _Caudovirales_. 



### Mapping and Tree Construction <a name="map"></a>


### Phage to Bacteria Ratio Calculation <a name="ratio"></a>





### Contact <a name="contact"></a>

Kristopher Kieft, kieft@wisc.edu
Zhichao Zhou 
