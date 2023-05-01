"""
This script was written by Kelsey Witt Dillon in July 2021. The goal of this script is to print a
list of archaic allele frequency differences between Han Chinese (CHB) individuals and individuals 
from a selected admixed American population (PEL, MXL, CLM, or PUR) from the 1000 genome project.
The script first identifies archaic alleles that are found in African populations at a frequency of 
less than 1% and a minimum frequency of 1% in the chosen American population, then calculates
the frequency difference for that allele between CHB and the AMR population.

The input file is gzipped csv file that lists the SNPs that are rare in Africa (<1% frequency) and 
shared with archaic individuals, and the genotypes for each individual: "archaic_rareAFR_genotypes.csv.gz",
These files are generated using the script "ID_rare_AFR_archaic_genotypes.py", which is accessed from a 
different github repository, github.com/kelsey-witt/archaic_analysis. ("infile", line 36)

Usage: python3 chb_amr_allele_freq_diff.py [AMR], where AMR is an admixed American population.
Possible inputs are "PEL", "CLM", "MXL", and "PUR".

This script generates two output files, with identical formats. The files are tab-delimited
text files with the following columns: Chromosome, Position, American allele frequency,
Han Chinese allele frequency, and the frequency difference (which is AMR frequency-CHB frequency).
The two outfiles are "archaic_snps_AMR_chbdiff_nr_1p.txt" and "archaic_snps_AMR_chbdiff_nr_1p.txt".
The _1p file includes all SNPs with a minimum AMR frequency of 1% while the _5p file includes
all SNPs with a minimum AMR frequency of 5%.
"""

import gzip
import sys

amrPop = sys.argv[1]

popTracker = {}
populations = [amrPop,"CHB"]

for pop in populations:
    popTracker[pop] = []

infile = "archaic_rareAFR_genotypes.csv.gz" 
outfile1 = "archaic_snps_" + amrPop + "_chbdiff_nr_1p.txt"
outfile5 = "archaic_snps_" + amrPop + "_chbdiff_nr_5p.txt"

def calc_total_X_allele(allele_total, genotype):
    allele_count_x = sum(genotype.count(x) for x in ("0", "1"))
    allele_total += allele_count_x
    return allele_total

def calc_pop_freq(population, nonafr_allele):
    pop_allele_count = 0
    allele_total = len(popTracker[population]) * 2
    if chr == "X":
        allele_total = 0
    for ind in popTracker[population]:
        pop_allele_count += spline[ind].count(nonafr_allele)
        if chrom == "X": #count number of X alleles as we go to account for recombinant and non-recombinant
            allele_total = calc_total_X_allele(allele_total, spline[ind])
    nonafr_frequency = float(pop_allele_count/allele_total)
    rounded_freq = round(nonafr_frequency, 3)
    return nonafr_frequency

g = open(outfile1, 'w')
h = open(outfile5, 'w')

outCols = ["chrom","pos","geno","alt","cha","vin","den","amr","chb","diff"]
outLine = "\t".join(outCols) + "\n"
g.write(outLine)
h.write(outLine)
with gzip.open(infile, 'rt') as f:
    for line in f:
        if "Chromosome" in line:
            spline = line[:-1].split(sep=",")
            for colnum in range(9,len(spline)):
                indPop = spline[colnum][0:3]
                if indPop in populations:
                    popTracker[indPop].append(colnum)
        else:
            spline = line[:-1].split(sep=",")
            chrom,pos,_,_,archAllele = spline[0:5]
            amrArchFreq = calc_pop_freq(amrPop,archAllele)
            if amrArchFreq > 0.01:
                chbArchFreq = calc_pop_freq("CHB",archAllele)
                freqDiff = amrArchFreq-chbArchFreq
                archaicHas = ["0","0","0","0"]
                archaicGenos = spline[5:9]
                for i in range(0,4):
                    if archAllele in archaicGenos[i]:
                        archaicHas[i] = "1"
                outCols = [str(chrom),str(pos),archAllele] + archaicHas + [str(amrArchFreq),str(chbArchFreq),str(freqDiff)]
                outLine = "\t".join(outCols) + "\n"
                g.write(outLine)
                if amrArchFreq > 0.05:
                    h.write(outLine)
