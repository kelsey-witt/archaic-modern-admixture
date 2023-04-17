"""
This script was written by Kelsey Witt Dillon in August 2018. The goal of this script is to calculate the number of 
archaic SNPs per ancestry region in admixed American populations. This script uses the 1000 Genomes data, including
the panel file, a series of bed files that identify the regions of European, Native American, and African ancestry
in the four admixed Native American populations, and the sprime calls for these individuals. The script first identifies
regions of the genome that have the same ancestry for both chromosomes, and then compares the list of archaic SNPs identified
by Sprime to the genotypes of the individuals in the 1000 Genomes dataset then tracks the number of archaic alleles
per individual and per ancestry type. This script can examine a variety of sets of archaic alleles, as defined by 
archaic_set, and six options are possible. "nd_either" includes all archaic alleles and "nd_both" includes archaic
alleles that are shared between Neanderthals and Denisovans (and identified as "match" in the Sprime data). "d_only"
and "n_only" refer to archaic SNPs that are unique to Denisovans or Neanderthals, respectively, and "d_all" and "n_all"
refer to all archaic SNPs that are found in Denisovans or Neanderthals, whether or not they are shared with the other
archaic human.
The output file is a csv file that lists the population, the individual, and the length of tracts (ie total number of 
base pairs of a given ancestry), number of archaic alleles, and density (# alleles/length) of each ancestry type.
"""
import gzip
import re
import sys
import csv

archaic_set = sys.argv[1]

#pop_file = "integrated_call_samples_v3.20130502.ALL.panel"
outfile = "arch_allele_counts_per_ind_ancestry_sprime_all_" + archaic_set + ".csv"
pop_file = "/users/kwittdil/data/data/modern_genomes/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
admixPops = ["PEL", "PUR", "CLM", "MXL"]
Admix_Pop_Tracker = {}

def decode_split_line(line):
    line = str(line.decode('utf-8'))
    if '#' not in line:
        line = line.split()
    return line

for pop in admixPops:
    Admix_Pop_Tracker[pop] = []
    

col_counter = 9
with open(pop_file, 'r') as p:
    next(p)
    for line in p:
        line_col = line.split()
        pop = line_col[1]
        if pop in admixPops:
            Admix_Pop_Tracker[pop].append(col_counter)
        col_counter += 1

popInd = {}
admixedIndA = {}
admixedIndB = {}
admixedIndRegions = {}
admixedIndArchSnps = {}
ancestries = ["AFR", "EUR", "NAT"]
nonafr_ancestries = ["EUR", "NAT"]
all_ancestries = ["AFR", "EUR", "NAT", "NonAfr"]
diplo_ancestries = ["AFR_AFR", "AFR_EUR", "AFR_NAT", "EUR_EUR", "EUR_NAT", "NAT_NAT"]
for pop in admixPops:
    popInd[pop] = []
with open(pop_file, 'r') as p:
    next(p)
    for line in p:
        line_col = line.split()
        if line_col[1] in admixPops:
            ind = line_col[0]
            pop = line_col[1]
            popInd[pop].append(ind)
            admixedIndRegions[ind] = {}
            admixedIndA[ind] = []
            admixedIndB[ind] = []
            for anc in diplo_ancestries:
                admixedIndRegions[ind][anc] = []
                #admixedIndArchSnps[ind][anc] = 0

#print(popInd)
#Dump ancestry regions into dictionaries for A and B

for pop in popInd:
    for ind in popInd[pop]:
        for genomeSide in ["A", "B"]:
            #infile = "./Ancestry/"+ pop + "/" + ind + "_" + genomeSide + "_final.bed"
            infile = "./ancestry_files/"+ pop + "/" + ind + "_" + genomeSide + "_final.bed"
            with open(infile) as f:
                for line in f:
                    spline = line.split()
                    chrom, regionStart, regionEnd, ancestry = spline[0:4]
                    if ancestry in ancestries and chrom != "X":  
                        regionStart = int(regionStart)
                        regionEnd = int(regionEnd)
                        chrom = int(chrom)
                        if genomeSide == "A":
                                admixedIndA[ind].append([chrom,regionStart,regionEnd,ancestry])
                        else:
                            admixedIndB[ind].append([chrom,regionStart,regionEnd,ancestry])

for pop in popInd:
    for ind in popInd[pop]:
        Apos = 0
        Bpos = 0
        while Apos < len(admixedIndA[ind]) and Bpos < len(admixedIndB[ind]):
            Aregion = admixedIndA[ind][Apos]
            Achrom, Astart, Aend, Aanc = Aregion[0:4]
            Bregion = admixedIndB[ind][Bpos]
            Bchrom, Bstart, Bend, Banc = Bregion[0:4]
            if Achrom == Bchrom:
                if Astart < Bend and Bstart < Aend:
                    regchrom = Achrom
                    regStart = max(Astart, Bstart)
                    regEnd = min(Aend, Bend)
                    ancList = [Aanc, Banc]
                    ancList.sort()
                    regAnc = ancList[0] + "_" + ancList[1]
                    admixedIndRegions[ind][regAnc].append([regchrom,regStart,regEnd,0])
                    if regEnd == Aend:
                        Apos += 1
                    else:
                        Bpos += 1
                elif Astart > Bend:
                    Bpos += 1
                elif Bstart > Aend:
                    Apos += 1
            elif Achrom < Bchrom:
                Apos +=1
            else:
                Bpos +=1


#Use previous archaic allele counter script to identify archaic alleles

chr_list = []
for c in range(1,23):
    chr_list.append(str(c))

for chromosome in chr_list:
    #modern_file = "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    modern_file = "/users/kwittdil/data/data/modern_genomes/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr" + chromosome + ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    alleleMatch = {}
    for pop in admixPops:
        sprime_file = "./sprime_data/" + pop + ".chr" + chromosome + ".ND_match"
        #sprime_file = "./Sprime/perchrom_files/" + pop + ".chr" + chromosome + ".ND_match"
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
                    if position not in alleleMatch.keys():
                        alleleMatch[position] = arch_allele
    with gzip.open(modern_file, 'r') as mf:
        for mod_line in mf:
            mod_line = decode_split_line(mod_line)
            mod_position = mod_line[1]
            if mod_position in alleleMatch.keys():
                for pop in Admix_Pop_Tracker:
                    indCounter = 0
                    arch_allele = alleleMatch[mod_position]
                    for ind in Admix_Pop_Tracker[pop]:
                        if arch_allele in mod_line[ind]:
                            individual = popInd[pop][indCounter]
                            for anc in diplo_ancestries:
                                #print(admixedIndRegions[individual])
                                for region in admixedIndRegions[individual][anc]:
                                    if str(region[0])==chromosome:
                                    #if the SNP is in a region of a given ancestry for that individual
                                        if int(mod_position) >= region[1] and int(mod_position) <= region[2]:
                                            region[3] += 1
                        indCounter += 1

with open(outfile, 'w', newline = '') as csvfile:
    w = csv.writer(csvfile, delimiter=",")
    header = ["pop_ancestry", "ind", "segment", "alleles", "len", "density"]
    w.writerow(header)
    for pop in admixPops:
        for ind in popInd[pop]:
            for anc in diplo_ancestries:
                ancestry_len = 0
                #print(str(len(admixedIndRegions[ind][anc])))
                if len(admixedIndRegions[ind][anc]) != 0:
                    for anc_region in admixedIndRegions[ind][anc]:
                        print(anc_region)
                        alleleCt = anc_region[3]
                        allele_line = []
                        allele_line.append(pop + "_" + anc)
                        allele_line.append(ind)
                        segment = anc_region[2] - anc_region[1]
                        coord = str(anc_region[0]) + ":" + str(anc_region[1]) + "-" + str(anc_region[2])
                        allele_density = alleleCt/segment
                        allele_line.append(coord)
                        allele_line.append(str(alleleCt))
                        allele_line.append(str(segment))
                        allele_line.append(str(allele_density))
                        w.writerow(allele_line)