"""
This script was written by Kelsey Witt Dillon in January 2019. The goal of this script is to calculate the number of 
archaic SNPs per individual in the 1000 Genomes Dataset, using the Sprime archaic calls to define archaic SNPs. 
This script uses the 1000 Genomes data, including the panel file, and the sprime calls for these individuals. The 
script first identifies which archaic alleles have been called using Sprime for the population, and then counts the 
total number of archaic SNPs for each individual. This script can examine a variety of sets of archaic alleles, as defined by 
archaic_set, and six options are possible. "nd_either" includes all archaic alleles and "nd_both" includes archaic
alleles that are shared between Neanderthals and Denisovans (and identified as "match" in the Sprime data). "d_only"
and "n_only" refer to archaic SNPs that are unique to Denisovans or Neanderthals, respectively, and "d_all" and "n_all"
refer to all archaic SNPs that are found in Denisovans or Neanderthals, whether or not they are shared with the other
archaic human.
The output file is a csv file that lists the individual, their population, the number of sites containing archaic SNPs 
and the total number of archaic SNPs across the genome.
"""

import gzip
import re
import sys
import csv

archaic_set = sys.argv[1] #possible options: d_only, n_only, d_all, n_all, nd_both, nd_either

#pop_file = "integrated_call_samples_v3.20130502.ALL.panel"
pop_file = "/home/1000genome_data/integrated_call_samples_v3.20130502.ALL.panel"
outfile = archaic_set + "_nonafr_allele_pos_ind_sprime_popspec_1902.csv"

populations = ["CHB", "JPT", "CHS", "CDX", "KHV", "CEU", "TSI", "FIN", "GBR", "IBS", "MXL", "PUR", "CLM", "PEL", "GIH", "PJL", "BEB", "STU", "ITU"]
Ind_Tracker = {}
Pop_Tracker = {}

for pop in populations:
    Pop_Tracker[pop] = {}

#chr_list = ["22"]
arch_allele_total = 0

chr_list = []
for c in range(1,23):
    chr_list.append(str(c))

col_counter = 9
with open(pop_file, 'r') as p:
    next(p)
    for line in p:
        line_col = line.split()
        pop = line_col[1]
        if pop in populations:
            Ind_Tracker[col_counter]=line_col[0]
            Pop_Tracker[pop][col_counter] = [0,0] #first number is positions, second number is alleles
        col_counter += 1

def decode_split_line(line):
    line = str(line.decode('utf-8'))
    if '#' not in line:
        line = line.split()
    return line

for chromosome in chr_list:
    #modern_file = "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    modern_file = "/home/1000genome_data/ALL.chr" + chromosome + ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    alleleMatch = {}
    allArchaicSNPs = set()
    for pop in populations:
        sprime_file = "./sprime_data/" + pop + ".chr" + chromosome + ".ND_match"
        #sprime_file = "./Sprime/perchrom_files/" + pop + ".chr" + chromosome + ".ND_match"
        alleleMatch[pop] = {}
        with open(sprime_file, 'r') as sf:
            next(sf)
            for line in sf:
                line_split = line.split()
                position = line_split[1]
                arch_allele = line_split[6] 
                Nstatus, Dstatus = line_split[8:10]
                isArchaic = False
                if Nstatus == "match" or Dstatus == "match":
                    if archaic_set == "nd_either":
                        isArchaic = True
                    elif archaic_set == "nd_both":
                        if Nstatus == "match" and Dstatus == "match":
                            isArchaic = True
                    elif archaic_set == "n_all":
                        if Nstatus == "match":
                            isArchaic = True
                    elif archaic_set == "n_only":
                        if Nstatus == "match" and Dstatus == "mismatch":
                            isArchaic = True
                    elif archaic_set == "d_all":
                        if Dstatus == "match":
                            isArchaic = True
                    elif archaic_set == "d_only":
                        if Dstatus == "match" and Nstatus == "mismatch":
                            isArchaic = True
                if isArchaic:
                    allArchaicSNPs.add(position)
                    alleleMatch[pop][position] = arch_allele

    with gzip.open(modern_file, 'r') as mf:
        for mod_line in mf:
            mod_line = decode_split_line(mod_line)
            mod_position = mod_line[1]
            if mod_position in allArchaicSNPs:
                #print(str(mod_position))
                for pop in populations:
                    if mod_position in alleleMatch[pop].keys():
                        arch_allele = alleleMatch[pop][mod_position]
                        for ind in Pop_Tracker[pop].keys():
                            if arch_allele in mod_line[ind]:
                                Pop_Tracker[pop][ind][0] += 1
                                num_archaic = mod_line[ind].count(arch_allele)
                                Pop_Tracker[pop][ind][1] += num_archaic

with open(outfile, 'w', newline = '') as csvfile:
    w = csv.writer(csvfile, delimiter=",")
    header = ["ind", "pop", "positions", "alleles"]
    w.writerow(header)
    for pop in populations:
        for ind in Pop_Tracker[pop].keys():
            indID = Ind_Tracker[ind]
            indPop = pop
            indPos = str(Pop_Tracker[pop][ind][0])
            indAllele = str(Pop_Tracker[pop][ind][1])
            score_line = [indID,indPop,indPos,indAllele]
            w.writerow(score_line)
