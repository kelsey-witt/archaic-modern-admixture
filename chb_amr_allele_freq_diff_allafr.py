"""
This script was written by Kelsey Witt Dillon in July 2021. The goal of this script is to print a
list of allele frequency differences between Han Chinese (CHB) individuals and individuals from
a selected admixed American population (PEL, MXL, CLM, or PUR) from the 1000 genome project.
The script first identifies alleles that are found in African populations at a frequency of 
less than 1% and a minimum frequency of 1% in the chosen American population, then calculates
the frequency difference for that allele between CHB and the AMR population.

The input file is a series of gzipped csv files that lists the SNPs that are rare in Africa 
(<1% frequency) and the genotypes for each individual: "rare_AFR_genotypes_chrN.csv.gz",
where N is the chromosome number ("infile", line 70). These files are generated using the script 
"ID_rare_AFR_snps_genos.py", which is accessed from a different github repository, 
github.com/kelsey-witt/archaic_analysis.

Usage: python3 chb_amr_allele_freq_diff_allafr.py [AMR], where AMR is an admixed American population.
Possible inputs are "PEL", "CLM", "MXL", and "PUR".

This script generates two output files, with identical formats. The files are tab-delimited
text files with the following columns: Chromosome, Position, American allele frequency,
Han Chinese allele frequency, and the frequency difference (which is AMR frequency-CHB frequency).
The two outfiles are "nonafr_snps_AMR_chbdiff_nr_1p.txt" and "nonafr_snps_AMR_chbdiff_nr_1p.txt".
The _1p file includes all SNPs with a minimum AMR frequency of 1% while the _5p file includes
all SNPs with a minimum AMR frequency of 5%.
"""

import gzip
import sys

amrPop = sys.argv[1]

popTracker = {}
populations = [amrPop,"CHB"]

chr_list = ["X"]
for c in range(1,23):
    chr_list.append(str(c))

for pop in populations:
    popTracker[pop] = []

outfile1 = "nonafr_snps_" + amrPop + "_chbdiff_nr_1p.txt"
outfile5 = "nonafr_snps_" + amrPop + "_chbdiff_nr_5p.txt"

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
        if chr == "X": #count number of X alleles as we go to account for recombinant and non-recombinant
            allele_total = calc_total_X_allele(allele_total, spline[ind])
    nonafr_frequency = float(pop_allele_count/allele_total)
    #rounded_freq = round(nonafr_frequency, 3)
    return nonafr_frequency

g = open(outfile1, 'w')
h = open(outfile5, 'w')

outCols = ["chrom","pos","amr","chb","diff"]
outLine = "\t".join(outCols) + "\n"
g.write(outLine)
h.write(outLine)
for chr in chr_list:
    infile = "rare_AFR_genotypes_chr" + chr + ".csv.gz"
    with gzip.open(infile, 'rt') as f:
        for line in f:
            if "Chromosome" in line:
                if chr == 'X':
                    spline = line[:-1].split(sep=",")
                    for colnum in range(9,len(spline)):
                        indPop = spline[colnum][0:3]
                        if indPop in populations:
                            popTracker[indPop].append(colnum)
            else:
                spline = line[:-1].split(sep=",")
                pos = spline[1]
                nonAfrAllele = spline[4]
                amrArchFreq = calc_pop_freq(amrPop,nonAfrAllele)
                if amrArchFreq > 0.01:
                    chbArchFreq = calc_pop_freq("CHB",nonAfrAllele)
                    freqDiff = amrArchFreq-chbArchFreq
                    outCols = [str(chr),str(pos),str(amrArchFreq),str(chbArchFreq),str(freqDiff)]
                    outLine = "\t".join(outCols) + "\n"
                    g.write(outLine)
                    if amrArchFreq > 0.05:
                        h.write(outLine)
