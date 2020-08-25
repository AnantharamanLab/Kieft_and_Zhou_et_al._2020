#! /usr/bin/env python3
# Author: Kristopher Kieft, UW-Madison, 2020

###############
## This script is associated with the manuscript
## Ecology of inorganic sulfur auxiliary metabolism in widespread bacteriophages
## Kieft and Zhou et al. 2020
###############

# Usage: $ python3 Kieft_and_Zhou_et_al_2020.phage-taxonomy.py -i <input_fasta_file> -t <threads> -f <format>
# PTT: phage taxonomy tool
# Version comment: the database was compiled using Diamond v2.0.0.138

# Contact Kristopher Kieft (kieft@wisc.edu) with questions regarding this script

import warnings
warnings.filterwarnings("ignore")
import sys
import argparse
import subprocess
import pandas as pd
from collections import Counter
from Bio.SeqIO.FastaIO import SimpleFastaParser
import os
import time
import logging
from datetime import date
import datetime

start = time.time()
round((time.time() - start)/60,1)
start_time = str(datetime.datetime.now().time()).rsplit(".",1)[0]

############################### Set Arguments  #################################
PTT_path = str(os.path.dirname(os.path.abspath(__file__)))
PTT = argparse.ArgumentParser(description='Estimates phage taxonomy.')
PTT.add_argument('--version', action='version', version='PTT v0.1.0')

PTT.add_argument('-i', type=str, nargs=1, required=True, help='input protein fasta file. Must be in format "name_#" with no spaces in protein name (Prodigal format).')
PTT.add_argument('-t', type=str, nargs=1, default='1', help='number of threads [default=1]')
PTT.add_argument('-f', type=str, nargs=1, default='prot', choices=['prot','nucl'], help='format of input [default="prot"]')
PTT.add_argument('-d', type=str, nargs=1, default=str(PTT_path) + '/PTT_database.dmnd', help='path to diamond database')
PTT.add_argument('-v', type=str, nargs=1, default=str(PTT_path) + '/PTT_virus_taxonomy.tsv', help='path to reference taxonomy file')

args = PTT.parse_args()
input = str(args.i[0])
database = args.d
taxonomy = args.v
threads = args.t[0]
format = str(args.f[0])
try:
    filepath = str(args.i[0].rsplit('/',1)[1])
    infile = str((filepath.rsplit('.',1)[:-1])[0])
except Exception:
    infile = str((input.rsplit('.',1)[:-1])[0])
base = infile.rsplit(".",1)[0]
############################### Set input  #####################################

if format == "prot" or format == "p":
    sequences = []
    with open(str(input), 'r') as fasta:
        with open(str(base)+'.temp.fasta', 'w') as rename:
            for name, seq in SimpleFastaParser(fasta):
                try:
                    temp = name.split(' # ',1)[0]
                    check = name.split(' # ',1)[1]
                    rename.write('>' + str(temp).replace(' ','~~') + '\n' + str(seq) + '\n')
                    sequences.append(str(temp).rsplit('_',1)[0])
                except Exception:
                    try:
                        temp = name.split('\t',1)[0]
                        check = name.split('\t',1)[1]
                        rename.write('>' + str(temp).replace(' ','~~') + '\n' + str(seq) + '\n')
                        sequences.append(str(temp).rsplit('_',1)[0])
                    except Exception:
                        rename.write('>' + str(name).replace(' ','~~') + '\n' + str(seq) + '\n')
                        sequences.append(str(name).rsplit('_',1)[0])
    sequences = list(set(sequences))
    subprocess.run("diamond blastp -d " + str(database) + " -q " + str(base)+'.temp.fasta' + " -o " + str(base) + ".PTT.diamond.out --max-target-seqs 3 -e 1e-5 --query-cover 10 --subject-cover 10 --threads " + str(threads) + " --quiet -f 6 qseqid sseqid pident evalue bitscore", shell=True)

if format == "nucl" or format == "n":
    sequences = []
    with open(str(input), 'r') as fasta:
        with open(str(base)+'.temp.fasta', 'w') as rename:
            for name, seq in SimpleFastaParser(fasta):
                rename.write('>' + str(name).replace(' ','~~') + '\n' + str(seq) + '\n')
                sequences.append(str(name))
    sequences = list(set(sequences))
    subprocess.run(['prodigal', '-i', str(input)+'.temp.fasta', '-a', str(base)+'.faa', '-p', 'meta', '-q', '-o', infile+'.prodigal.temp'])
    with open(str(base)+'.faa', 'r') as faa:
        with open(str(base)+'.temp.faa', 'w') as temp_faa:
            for name, seq in SimpleFastaParser(faa):
                temp = name.split(' # ',1)[0]
                temp_faa.write('>' + str(temp) + '\n' + str(seq) + '\n')
    subprocess.run("diamond blastp -d " + str(database) + " -q " + str(base)+'.temp.faa' + " -o " + str(base) + ".PTT.diamond.out --max-target-seqs 3 -e 1e-5 --query-cover 10 --subject-cover 10 --threads " + str(threads) + " --quiet -f 6 qseqid sseqid pident evalue bitscore", shell=True)

subprocess.run('cat ' + str(base) + '.PTT.diamond.out | sed "s/~~/ /g" 1> ' + str(base) + '.PTT.diamond 2>/dev/null', shell=True)
subprocess.run('rm ' + str(base) + '.PTT.diamond.out 2>/dev/null', shell=True)
subprocess.run('rm ' + str(base)+'.temp.fasta 2>/dev/null', shell=True)
subprocess.run('rm ' + str(base)+'.temp.faa 2>/dev/null', shell=True)
subprocess.run('rm ' + str(base)+'.prodigal.temp 2>/dev/null', shell=True)

with open(str(taxonomy), "r") as reffile:
    n = 5
    ref_dict = {}
    ref_dict_genus = {}
    ref = reffile.read().replace("\n","\t").split("\t")
    while n < len(ref):
        ref_dict.update({ref[n]:str(ref[n+1])+'~'+str(ref[n+2])+'~'+str(ref[n+3])})
        ref_dict_genus.update({ref[n]:str(ref[n+1])+'~'+str(ref[n+2])+'~'+str(ref[n+3])+'~'+str(ref[n+4])})
        n += 5

with open(str(base) + ".PTT.diamond", "r") as diamond:
    with open(str(base) + ".PTT.virus-topreps.tsv", "w") as outfile:
        with open(str(base) + ".PTT.virus-tophits.tsv", "w") as temp:
            diamond = diamond.read().replace("\n","\t").split("\t")
            if diamond[-1] == '':
                diamond = diamond[:-1]
            diamond.extend(['\t','\t','\t','\t','\t'])
            outfile.write("scaffold\tconf_1\to_1\tf_1\tsf_1\tconf_2\to_2\tf_2\tsf_2\tconf_3\to_3\tf_3\tsf_3\n")
            temp.write("scaffold\to_1\to_1v\to_2\to_2v\tf_1\tf_1v\tf_2\tf_2v\tsf_1\tsf_1v\tsf_2\tsf_2v\n")
            n = 0
            new = True
            tax_est = []
            order_list = []
            family_list = []
            subfamily_list = []
            while n < len(diamond)-5:
                if new == True:
                    outfile.write(str(diamond[n].rsplit("_",1)[0])+"\t")
                    temp.write(str(diamond[n].rsplit("_",1)[0])+"\t")
                    new = False
                if float(diamond[n+2]) >= 0.2:
                    prot = diamond[n+1].rsplit("_",1)[0]
                    order_list.append(ref_dict[prot].split("~")[0])
                    family_list.append(ref_dict[prot].split("~")[1])
                    subfamily_list.append(ref_dict[prot].split("~")[2])
                    tax_est.append(ref_dict[prot])
                if str(diamond[n].rsplit("_",1)[0]) != str(diamond[n+5].rsplit("_",1)[0]):
                    order_count = Counter(order_list)
                    family_count = Counter(family_list)
                    subfamily_count = Counter(subfamily_list)
                    order1 = str(order_count.most_common(1)[0]).replace(")","").replace("(","").replace("'","").split(",")[0]
                    order1v = str(order_count.most_common(1)[0]).replace(")","").replace("(","").replace("'","").split(",")[1]
                    family1 = str(family_count.most_common(1)[0]).replace(")","").replace("(","").replace("'","").split(",")[0]
                    family1v = str(family_count.most_common(1)[0]).replace(")","").replace("(","").replace("'","").split(",")[1]
                    subfamily1 = str(subfamily_count.most_common(1)[0]).replace(")","").replace("(","").replace("'","").split(",")[0]
                    subfamily1v = str(subfamily_count.most_common(1)[0]).replace(")","").replace("(","").replace("'","").split(",")[1]
                    try:
                        order2 = str(order_count.most_common()[1]).replace(")","").replace("(","").replace("'","").split(",")[0]
                        order2v = str(order_count.most_common()[1]).replace(")","").replace("(","").replace("'","").split(",")[1]
                    except Exception:
                        order2 = 'unassigned'
                        order2v = '0'
                    try:
                        family2 = str(family_count.most_common()[1]).replace(")","").replace("(","").replace("'","").split(",")[0]
                        family2v = str(family_count.most_common()[1]).replace(")","").replace("(","").replace("'","").split(",")[1]
                    except Exception:
                        family2 = 'unassigned'
                        family2v = '0'
                    try:
                        subfamily2 = str(subfamily_count.most_common()[1]).replace(")","").replace("(","").replace("'","").split(",")[0]
                        subfamily2v = str(subfamily_count.most_common()[1]).replace(")","").replace("(","").replace("'","").split(",")[1]
                    except Exception:
                        subfamily2 = 'unassigned'
                        subfamily2v = '0'
                    temp.write(str(order1) + "\t" + str(order1v) + "\t" + str(order2) + "\t" + str(order2v) + "\t" + str(family1) + "\t" + str(family1v) + "\t" + str(family2) + "\t" + str(family2v) + "\t" + str(subfamily1) + "\t" + str(subfamily1v) + "\t" + str(subfamily2) + "\t" + str(subfamily2v) + "\n")
                    order_list = []
                    family_list = []
                    subfamily_list = []
                    genus_list = []
                    tax_count = Counter(tax_est)
                    tax1 = str(tax_count.most_common(1)[0]).replace(")","").replace("(","").replace("'","").split(",")[0]
                    tax1val = str(tax_count.most_common(1)[0]).replace(")","").replace("(","").replace("'","").split(",")[1]
                    try:
                        tax2 = str(tax_count.most_common()[1]).replace(")","").replace("(","").replace("'","").split(",")[0]
                        tax2val = str(tax_count.most_common()[1]).replace(")","").replace("(","").replace("'","").split(",")[1]
                    except Exception:
                        tax2 = 'unassigned~unassigned~unassigned'
                        tax2val = '0'
                    try:
                        tax3 = str(tax_count.most_common()[2]).replace(")","").replace("(","").replace("'","").split(",")[0]
                        tax3val = str(tax_count.most_common()[2]).replace(")","").replace("(","").replace("'","").split(",")[1]
                    except Exception:
                        tax3 = 'unassigned~unassigned~unassigned'
                        tax3val = '0'
                    outfile.write(str(tax1val) + "\t" + str(tax1).split("~")[0] + "\t" + str(tax1).split("~")[1] + "\t" + str(tax1).split("~")[2] + "\t" + str(tax2val) + "\t" + str(tax2).split("~")[0] + "\t" + str(tax2).split("~")[1] + "\t" + str(tax2).split("~")[2] + "\t" + str(tax3val) + "\t" + str(tax3).split("~")[0] + "\t" + str(tax3).split("~")[1] + "\t" + str(tax3).split("~")[2] + "\n")
                    new = True
                    tax_est = []
                n += 5

written = []
with open(str(base) + ".PTT.virus-topreps.tsv", "r") as taxfile:
    with open(str(base) + ".PTT.virus-taxonomy.tsv", "w") as outfile:
        with open(str(base) + ".PTT.virus-tophits.tsv", "r") as helper:
            taxfile = taxfile.read().replace("\n","\t").split("\t")
            helper = helper.read().split("\n")
            if helper[-1] == '':
                helper = helper[:-1]
            helper_dict = {}
            for item in helper:
                h1 = item.split("\t",1)[0]
                h2 = item.split("\t",1)[1]
                helper_dict.update({str(h1):str(h2).replace("\t","~")})
            if taxfile[-1] == '':
                taxfile = taxfile[:-1]
            outfile.write("scaffold\tpred_score\torder\tfamily\tsubfamily\n")
            n = 13
            while n < len(taxfile):
                order = 'ambiguous'
                family = 'ambiguous'
                subfamily = 'ambiguous'
                oscore = 1
                fscore = 1
                sfscore = 1
                i = 0
                if float(taxfile[n+1]) >= 3:
                    if float(taxfile[n+1]) > float(taxfile[n+5]):
                        if str(helper_dict[str(taxfile[n])].split("~")[0]) == str(taxfile[n+2]):
                            order = str(taxfile[n+2])
                        elif str(taxfile[n+6]) == str(taxfile[n+2]):
                            order = str(taxfile[n+2])
                        else:
                            order = 'ambiguous'
                        if str(helper_dict[str(taxfile[n])].split("~")[4]) == str(taxfile[n+3]) and order != 'ambiguous' and str(order) == str(taxfile[n+2+i]):
                            family = str(taxfile[n+3])
                        elif str(taxfile[n+7]) == str(taxfile[n+3]) and order != 'ambiguous' and str(order) == str(taxfile[n+2+i]):
                            family = str(taxfile[n+3])
                        else:
                            family = 'ambiguous'
                        if str(helper_dict[str(taxfile[n])].split("~")[8]) == str(taxfile[n+4]) and family != 'ambiguous' and str(family) == str(taxfile[n+3+i]):
                            subfamily = str(taxfile[n+4])
                        elif str(taxfile[n+8]) == str(taxfile[n+4]) and family != 'ambiguous' and str(family) == str(taxfile[n+3+i]):
                            subfamily = str(taxfile[n+4])
                        else:
                            subfamily = 'ambiguous'
                    if float(taxfile[n+1]) == float(taxfile[n+5]):
                        if str(helper_dict[str(taxfile[n])].split("~")[0]) == str(taxfile[n+2]):
                            order = str(taxfile[n+2])
                        elif str(helper_dict[str(taxfile[n])].split("~")[0]) == str(taxfile[n+6]):
                            order = str(taxfile[n+6])
                            i = 4
                        elif str(taxfile[n+6]) == str(taxfile[n+2]):
                            order = str(taxfile[n+2])
                        else:
                            order = 'ambiguous'
                        if str(helper_dict[str(taxfile[n])].split("~")[4]) == str(taxfile[n+3+i]) and order != 'ambiguous':
                            if float(helper_dict[str(taxfile[n])].split("~")[5]) == float(helper_dict[str(taxfile[n])].split("~")[7]):
                                if str(taxfile[n+7+i]) == 'unassigned' and str(order) == str(taxfile[n+2+i]):
                                    family = str(taxfile[n+3+i])
                                elif str(taxfile[n+3+i]) == 'unassigned' and str(order) == str(taxfile[n+6+i]):
                                    family = str(taxfile[n+7+i])
                                else:
                                    family = 'ambiguous'
                            elif float(helper_dict[str(taxfile[n])].split("~")[5]) > float(helper_dict[str(taxfile[n])].split("~")[7]):
                                if str(helper_dict[str(taxfile[n])].split("~")[4]) == str(taxfile[n+3+i]) and str(order) == str(taxfile[n+2+i]):
                                    family = str(taxfile[n+3+i])
                                elif str(helper_dict[str(taxfile[n])].split("~")[4]) == str(taxfile[n+7+i]) and str(order) == str(taxfile[n+6+i]):
                                    family = str(taxfile[n+7+i])
                                else:
                                    family = 'ambiguous'
                            elif str(order) == str(taxfile[n+2+i]):
                                family = str(taxfile[n+3+i])
                            else:
                                family = 'ambiguous'
                        elif str(helper_dict[str(taxfile[n])].split("~")[4]) == str(taxfile[n+7+i]) and order != 'ambiguous' and float(helper_dict[str(taxfile[n])].split("~")[5]) >= float(helper_dict[str(taxfile[n])].split("~")[7]):
                            if float(taxfile[n+1+i]) == float(taxfile[n+5+i]):
                                if str(taxfile[n+7+i]) == 'unassigned' and str(order) == str(taxfile[n+2+i]):
                                    family = str(taxfile[n+3+i])
                                elif str(taxfile[n+3+i]) == 'unassigned' and str(order) == str(taxfile[n+6+i]):
                                    family = str(taxfile[n+7+i])
                                else:
                                    family = 'ambiguous'
                            elif str(order) == str(taxfile[n+6+i]):
                                family = str(taxfile[n+7+i])
                            else:
                                family = 'ambiguous'
                        elif str(taxfile[n+7+i]) == str(taxfile[n+3]) and order != 'ambiguous' and str(order) == str(taxfile[n+2+i]):
                            family = str(taxfile[n+3+i])
                        elif str(taxfile[n+7]) == str(taxfile[n+11]) and order != 'ambiguous' and str(order) == str(taxfile[n+6+i]):
                            family = str(taxfile[n+7+i])
                        else:
                            family = 'ambiguous'
                        if str(helper_dict[str(taxfile[n])].split("~")[8]) == str(taxfile[n+4+i]) and family != 'ambiguous' and str(family) == str(taxfile[n+3+i]):
                            subfamily = str(taxfile[n+4+i])
                        elif str(taxfile[n+8-i]) == str(taxfile[n+4+i]) and family != 'ambiguous' and str(family) == str(taxfile[n+3+i]):
                            subfamily = str(taxfile[n+4+i])
                        elif str(helper_dict[str(taxfile[n])].split("~")[8]) == str(taxfile[n+8+i]) and family != 'ambiguous' and str(family) == str(taxfile[n+7+i]):
                            subfamily = str(taxfile[n+8+i])
                        else:
                            subfamily = 'ambiguous'
                    if str(taxfile[n+2]) == str(order):
                        oscore += 1
                    if str(taxfile[n+6]) == str(order):
                        oscore += 1
                    if str(taxfile[n+10]) == str(order):
                        oscore += 1
                    if str(taxfile[n+3]) == str(family):
                        fscore += 1
                    if str(taxfile[n+7]) == str(family):
                        fscore += 1
                    if str(taxfile[n+11]) == str(family):
                        fscore += 1
                    if str(taxfile[n+4]) == str(subfamily):
                        sfscore += 1
                    if str(taxfile[n+8]) == str(subfamily):
                        sfscore += 1
                    if str(taxfile[n+12]) == str(subfamily):
                        sfscore += 1
                    if str(helper_dict[str(taxfile[n])].split("~")[0]) == str(order):
                        oscore += 1
                    if str(helper_dict[str(taxfile[n])].split("~")[2]) == str(order):
                        oscore += 1
                    if str(helper_dict[str(taxfile[n])].split("~")[4]) == str(family):
                        fscore += 1
                    if str(helper_dict[str(taxfile[n])].split("~")[6]) == str(family):
                        fscore += 1
                    if str(helper_dict[str(taxfile[n])].split("~")[8]) == str(subfamily):
                        sfscore += 1
                    if str(helper_dict[str(taxfile[n])].split("~")[10]) == str(subfamily):
                        sfscore += 1
                    score = str(oscore)+"~"+str(fscore)+"~"+str(sfscore)
                else:
                    order = 'unknown'
                    family = 'unknown'
                    subfamily = 'unknown'
                    score = '0~0~0~0'
                written.append(str(taxfile[n]))
                outfile.write(str(taxfile[n]) + "\t" + str(score) + "\t" + str(order) + "\t" + str(family) + "\t" + str(subfamily) + "\n")
                n += 13

            for name in sequences:
                if name not in written:
                    outfile.write(str(name) + '\t0~0~0~0\tunknown\tunknown\tunknown\n')

with open(str(base) + ".PTT.diamond", "r") as diamond:
    table = pd.read_csv(diamond, sep="\t", header=None, names=["query", "subject", "pident", "evalue", "score"])
    sort = table.sort_values(by='evalue', ascending=True)
    drop = sort.drop_duplicates(subset='query', keep='first')
    write = drop.to_csv(str(base) + ".PTT.diamond-parsed.tsv", index=False, sep="\t")

with open(str(base) + ".PTT.diamond-parsed.tsv", "r") as diamond:
    with open(str(base) + ".PTT.protein-taxonomy.tsv", "w") as outfile:
        diamond = diamond.read().replace("\n","\t").split("\t")
        if diamond[-1] == '':
            diamond = diamond[:-1]
        diamond.extend(['\t','\t','\t','\t','\t'])
        outfile.write("scaffold\torder\tfamily\tsubfamily\tgenus\n")
        n = 5
        while n < len(diamond)-5:
            outfile.write(str(diamond[n])+"\t")
            if float(diamond[n+2]) >= 0.2:
                prot = diamond[n+1].rsplit("_",1)[0]
                order = ref_dict_genus[prot].split("~")[0]
                family = ref_dict_genus[prot].split("~")[1]
                subfamily = ref_dict_genus[prot].split("~")[2]
                genus = ref_dict_genus[prot].split("~")[3]
                outfile.write(str(order) + "\t" + str(family) + "\t" + str(subfamily) + "\t" + str(genus) + "\n")
            n += 5

#subprocess.run("rm " + str(base) + ".PTT.virus-tophits.tsv", shell=True)

logging.basicConfig(filename='PTT_log_' + str(base) + '.log', level=logging.INFO, format='%(message)s')
log_command = sys.argv
logging.info("Command:  " + str(" ".join(log_command)))
logging.info("Date:     " + str(date.today()))
logging.info("Start:    " + str(start_time))
logging.info("End:      " + str(datetime.datetime.now().time()).rsplit(".",1)[0])
logging.info("Runtime:  " + str(round((time.time() - float(start))/60,1)) + " minutes")
logging.info("Program:  PTT v0.1.0")
logging.info('                                                               ##')
logging.info('                                                             ##  ##')
logging.info('                                                           ##      ##')
logging.info('######   ##  ##     ##     #######   ######    #####       ##      ##')
logging.info('##  ##   ##  ##   ##  ##   ##        ##       ##             ##  ##')
logging.info('######   ######   ######   ##  ###   ######    ###             ##')
logging.info('##       ##  ##   ##  ##   ##   ##   ##           ##           ##')
logging.info('##       ##  ##   ##  ##   #######   ######   #####            ##')
logging.info('                                                            #  ##  #')
logging.info('                                                           # # ## # #')
logging.info('                                                          #   #  #   #')
logging.info('                                                         #            #')
logging.info("\n")

subprocess.run('mkdir PTT_' + str(base) + '_' + str(date.today()), shell=True)
subprocess.run('mkdir PTT_' + str(base) + '_' + str(date.today()) + '/analysis_files', shell=True)
subprocess.run('mv PTT_log_' + str(base) + '.log PTT_' + str(base) + '_' + str(date.today()), shell=True)
subprocess.run('mv ' + str(base) + '.PTT.*-taxonomy.tsv PTT_' + str(base) + '_' + str(date.today()), shell=True)
subprocess.run('mv ' + str(base) + '.PTT.* PTT_' + str(base) + '_' + str(date.today()) + '/analysis_files', shell=True)
if format == "nucl" or format == "n":
    subprocess.run('mv ' + str(base) + '.faa PTT_' + str(base) + '_' + str(date.today()), shell=True)
#
#
#
