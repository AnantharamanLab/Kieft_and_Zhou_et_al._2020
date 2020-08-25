# Kieft_and_Zhou_et_al._2020

This GitHub repository contains supplementary files associated with the manuscript _Ecology of inorganic sulfur auxiliary metabolism in widespread bacteriophages_

The bioRxiv preprint of this manuscript can be found [here](https://www.biorxiv.org/content/10.1101/2020.08.24.253096v1). 

8/25/2020  
Zhichao Zhou  
Kristopher Kieft  
Karthik Anantharaman  
University of Wisconsin-Madison  



### Explanations

##### Taxonomy
This folder contains the custom Python program and associated database to estimate the taxonomy of prokaryotic viruses. The program (`PTT.py`, Phage Taxonomy Tool) takes either proteins or genomes as an input file and will run Diamond Blastp against a database of reference viruses from NCBI GenBank. Diamond Blastp hits will be filtered and assignment of a taxonomic lineage will be determined according to counting Diamond Blastp hits. Assignments will be made on the levels of Order, Family and Sub-family, though database metadata is provided to the Genus level if manual curation of Genus is desired. Each taxonomic level will be assigned individually, starting with Order. If an Order is assigned it will move to Family before assigning Sub-family. For example, if Order _Caudovirales_ and Family _unknown_ is assigned then there will be no assignment of Sub-family since the higher level (i.e., Family) was unknown.  

_Dependencies_: Diamond, Prodigal (if input is nucleotides), Blast, BioPython. 

_Before starting_: Please unzip the Diamond database file using the command `gunzip PTT_database.dmnd.gz`. In the event you are using an uncompatible version of Diamond (v2 used here) you may need to generate the database on your own. Please unzip the protein database using the command `gunzip PTT_database.faa.gz` followed by the command `diamond makedb --in PTT_database.faa --db PTT_database.dmnd`. In addition to the Diamond database there is also a taxonomy metadata file required for the program (`PTT_virus_taxonomy.tsv`). Please place the database files (`PTT_database.dmnd` and `PTT_virus_taxonomy.tsv`) in the same location as the program script (`PTT.py`). If the databases are not in the same location you may choose to use the optional flags `-d` and `-v` to designate the database location(s).  

_Usage_: This programs takes a minimum of one argument (the input file). For quick usage, run `python3 PTT.py -i input_file.faa` with `input_file.faa` being a protein fasta file of your choice in Prodigal format. `-t` to increase the number of threads used. `-f` to indicate "prot" (protein) or "nucl" (nucleotide) input. `-d` to provide the location for the Diamond database (if not in same folder as `PTT.py`). `-v` to provide the location for the metadata taxonomy file (if not in same folder as `PTT.py`). See `PTT.py -h` for the help menu.  

_Outputs_: 

_CAUTION_: The input file headers (sequence names) cannot have spaces in the names or else the output will cut the names at those spaces. Input proteins must be in Prodigal format (e.g., >sequence_name_#).  
