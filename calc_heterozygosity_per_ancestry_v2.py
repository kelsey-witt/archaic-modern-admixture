"""
This script was written by Kelsey Witt Dillon in August 2019. The goal of this script is to calculate the heterozygosity
 per ancestry region in admixed American populations. This script uses the 1000 Genomes data, including
the panel file, a series of bed files that identify the regions of European, Native American, and African ancestry
in the four admixed Native American populations. The script first identifies regions of the genome that have the same 
ancestry for both chromosomes, and then counts the number of heterozygotes and homozygotes per region.
The output file is a csv file that lists the population, the individual, the chromosome, the number of region, and the
counts of heterozygotes and homozygotes
"""
import gzip
import re
import sys
import csv

#archaic_set = sys.argv[1]

#pop_file = "integrated_call_samples_v3.20130502.ALL.panel"
outfile = "heterozygosity_by_ancestry.csv"
pop_file = "/users/kwittdil/data/data/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
admixPops = ["PEL", "PUR", "CLM", "MXL"]
#pelIND = ["HG01572","HG01578","HG01927","HG01935","HG01965","HG01970","HG01979","HG01991","HG01992","HG02298","HG02304"]
#admixPops = ["PUR"]
Admix_Pop_Tracker = {}

def decode_split_line(line):
    line = str(line.decode('utf-8'))
    if '#' not in line:
        line = line.split()
    return line

def snp_is_biallelic(line):
    ref_allele = line[3]
    alt_allele = line[4]
    if len(ref_allele) == 1 and len(alt_allele) == 1:
        return True
    else:
        return False

for pop in admixPops:
    Admix_Pop_Tracker[pop] = []

col_counter = 9
with open(pop_file, 'r') as p:
    next(p)
    for line in p:
        line_col = line.split()
        pop = line_col[1]
        ind = line_col[0]
        if pop in admixPops:
            Admix_Pop_Tracker[pop].append(col_counter)
        col_counter += 1

popInd = {}
admixedIndA = {}
admixedIndB = {}
admixedIndRegions = {}
admixedIndRegionCounter = {}
admixedIndHet = {}
nonunk = ["AFR", "EUR", "NAT"]
#nonunk = ["EUR"]
ancestries = ["AFR_AFR", "AFR_EUR", "AFR_NAT", "EUR_EUR", "EUR_NAT", "NAT_NAT"]
#ancestries = ["EUR_EUR"]
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
            admixedIndRegionCounter[ind] = {}
            admixedIndHet[ind] = {}
            for anc in ancestries:
                admixedIndRegions[ind][anc] = []
                admixedIndRegionCounter[ind][anc] = 0
                admixedIndHet[ind][anc] = 0

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
                    if ancestry in nonunk and chrom != "X":  
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
                    admixedIndRegions[ind][regAnc].append([regchrom,regStart,regEnd])
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

chr_list = []
for c in range(1,23):
    chr_list.append(str(c))

#chr_list = ["22"]

for chromosome in chr_list:
    #modern_file = "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    modern_file = "/users/kwittdil/data/data/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr" + chromosome + ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    with gzip.open(modern_file, 'r') as mf:
        for mod_line in mf:
            mod_line = decode_split_line(mod_line)
            if '#' not in mod_line and snp_is_biallelic(mod_line):
                mod_position = mod_line[1]
                #print(mod_line)
                for pop in Admix_Pop_Tracker:
                    indCounter = 0
                    for ind in Admix_Pop_Tracker[pop]:
                        individual = popInd[pop][indCounter]
                        for anc in ancestries:
                            if len(admixedIndRegions[individual][anc])>0 and admixedIndRegionCounter[individual][anc]<len(admixedIndRegions[individual][anc]):
                                regionChr, regionSt, regionEnd = admixedIndRegions[individual][anc][admixedIndRegionCounter[individual][anc]][0:3]
                                while int(regionChr) < int(chromosome) and admixedIndRegionCounter[individual][anc]<(len(admixedIndRegions[individual][anc]))-1:
                                    admixedIndRegionCounter[individual][anc] += 1
                                    regionChr, regionSt, regionEnd = admixedIndRegions[individual][anc][admixedIndRegionCounter[individual][anc]][0:3]
                                if str(regionChr)==str(chromosome):
                                    #print("inchr")
                                    inChr=True
                            #if the SNP is in a region of a given ancestry for that individual
                                    if int(mod_position) >= regionSt and int(mod_position) < regionEnd:
                                        inRegion=True
                                        inChr = True
                                        if mod_line[ind] == "0|1" or mod_line[ind] == "1|0":
                                            #print("het")
                                            admixedIndHet[individual][anc] += 1
                                    if int(mod_position) >= int(regionEnd) and inRegion:
                                        admixedIndRegions[individual][anc][admixedIndRegionCounter[individual][anc]].append(admixedIndHet[individual][anc])
                                        admixedIndRegionCounter[individual][anc] += 1
                                        admixedIndHet[individual][anc] = 0
                                        inRegion=False
                                    elif int(mod_position) >= int(regionEnd) and inChr:
                                        admixedIndRegions[individual][anc][admixedIndRegionCounter[individual][anc]].append(admixedIndHet[individual][anc])
                                        admixedIndRegionCounter[individual][anc] += 1
                                        admixedIndHet[individual][anc] = 0
                                        inRegion=False 
                                        inChr = False
                        indCounter += 1
                #print(str(mod_position))
for pop in Admix_Pop_Tracker:
    indCounter = 0
    for ind in Admix_Pop_Tracker[pop]:
        individual = popInd[pop][indCounter]
        for anc in ancestries:
            if admixedIndRegionCounter[individual][anc]<len(admixedIndRegions[individual][anc]):
                admixedIndRegions[individual][anc][admixedIndRegionCounter[individual][anc]].append(admixedIndHet[individual][anc])
        indCounter += 1
                                
#print(admixedIndRegions)

with open(outfile, 'w', newline = '') as csvfile:
    w = csv.writer(csvfile, delimiter=",")
    header = ["pop", "ancestry", "ind", "region_chr", "num", "region_start", "region_end", "heterozygotes"]
    w.writerow(header)
    for pop in admixPops:
        for ind in popInd[pop]:
            for anc in ancestries:
                if len(admixedIndRegions[ind][anc])>0:
                    posCounter = 0
                    for region in admixedIndRegions[ind][anc]:
                        #print(region)
                        allele_line = [pop, anc, ind]
                        pulledInfo = region
                        regionNum = posCounter + 1
                        if len(region) > 0:
                            rChr, rSt, rEnd, rHet = pulledInfo[0:5]
                        else:
                            rChr = rSt = rEnd = rHet = "NA"
                        allele_line.append(str(rChr))
                        allele_line.append(str(regionNum))
                        allele_line.append(str(rSt))
                        allele_line.append(str(rEnd))
                        allele_line.append(str(rHet))
                        w.writerow(allele_line)
                        posCounter += 1
