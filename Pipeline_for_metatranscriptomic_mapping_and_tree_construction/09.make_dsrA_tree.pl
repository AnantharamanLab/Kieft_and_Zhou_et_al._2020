#!/usr/bin/perl

use strict;
use warnings;


#`cat dsrA_total.faa bacterial_and_viral_top_10_hit_accession_id.faa dsrA_ref.fasta > dsrA_total_and_ref.faa`;

`mafft dsrA_total_and_ref.faa > dsrA_total_and_ref.mafft.faa`;

`iqtree -nt AUTO -m MFP -bb 1000 -s dsrA_total_and_ref.mafft.faa -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -pre  dsrA_total_and_ref.mafft.faa -nt 60`;
