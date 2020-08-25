#!/usr/bin/perl

use strict;
use warnings;


`cat Lau_dsrA_protein.faa bacterial_and_viral_top_10_hit_accession_id.faa dsrA_ref.fasta > Lau_dsrA_protein_and_ref.faa`;

`mafft Lau_dsrA_protein_and_ref.faa > Lau_dsrA_protein_and_ref.mafft.faa`;

`iqtree -nt AUTO -m MFP -bb 1000 -s Lau_dsrA_protein_and_ref.mafft.faa -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -pre  Lau_dsrA_protein_and_ref.mafft.faa -nt 60`;
