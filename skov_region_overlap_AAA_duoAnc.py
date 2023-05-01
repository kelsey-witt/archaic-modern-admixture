"""
This script was written by Kelsey Witt Dillon in January 2023. The goal of this script is to 
identify how archaic introgression tracts from Skov et al. 2018 overlap with modern
ancestry tracts in admixed American individuals from Martin et al. 2017. The genomes of the 
individuals are subdivided into regions of diploid ancestry (ie AFR_AFR, EUR_AFR, etc.) and 
then the percentage of that ancestry tract that overlaps with a diploid archaic ancestry tract 
is calculated.

The input files are the .panel file from the 1000 genomes dataset ("pop_file", line 28),
the ancestry bed files from Martin et al., ("infile", line 76) and the list of diploid archaic
segments from Skov ("skovInfile", line 110)

Usage: python3 skov_region_overlap_AAA_duoAnc.py

This script generates a tab-delimited outfile, currently named "skov_ancestry_overlap.txt" ("outfile", 
line 66). The columns are the individual, the population, the diploid ancestry of the region, the 
chromosome, the start of the ancestry region, the end of the ancestry region, the length of the region
that includes archaic ancestry, and the percentage that length represents.
"""

import gzip
import sys

pop_file = "/users/kwittdil/data/data/modern_genomes/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
admixPops = ["PEL", "PUR", "CLM", "MXL"]
#chromosome = sys.argv[1]
"""
pop_file = "./integrated_call_samples_v3.20130502.ALL.panel"
admixPops = ["PUR"]
chromosome = "1"
"""

chromosome_list = []
for i in range(1,23):
    chromosome_list.append(str(i))

Admix_Pop_Tracker = {}

popInd = {}
admixedIndA = {}
admixedIndB = {}
admixedIndRegions = {}
admixedIndArchSnps = {}
ancestries = ["AFR", "EUR", "NAT"]
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

outputFile = "./skov_ancestry_overlap.txt"
f = open(outputFile, 'w') 
header = ["Individual", "Population", "Ancestry", "Chromosome", "Anc Start", "Anc End", "Shared Length", "% Shared"]
headerLine = "\t".join(header) + "\n"
f.write(headerLine)

for chromosome in chromosome_list:
#print(popInd)
#Sort through ancestry regions and pull list where both chromosomes have same ancestry
    for pop in popInd:
        for ind in popInd[pop]:
            for genomeSide in ["A", "B"]:
                #infile = "./" + pop + "/" + ind + "_" + genomeSide + "_final.bed"
                infile = "./ancestry_files/"+ pop + "/" + ind + "_" + genomeSide + "_final.bed"
                with open(infile) as j:
                    for line in j:
                        line_split = line.split()
                        chrom, regionStart, regionEnd, ancestry = line_split[0:4]
                        regionStart = int(regionStart)
                        regionEnd = int(regionEnd)
                        if str(chrom) == str(chromosome) and ancestry in ancestries:
                            if genomeSide == "A":
                                admixedIndA[ind].append([regionStart,regionEnd,ancestry])
                            else:
                                admixedIndB[ind].append([regionStart,regionEnd,ancestry])
            Apos = 0
            Bpos = 0
            while Apos < len(admixedIndA[ind]) and Bpos < len(admixedIndB[ind]):
                Aregion = admixedIndA[ind][Apos]
                Astart, Aend, Aanc = Aregion[0:4]
                Bregion = admixedIndB[ind][Bpos]
                Bstart, Bend, Banc = Bregion[0:4]
                if Astart < Bend and Bstart < Aend:
                    regStart = max(Astart, Bstart)
                    regEnd = min(Aend, Bend)
                    ancList = [Aanc, Banc]
                    ancList.sort()
                    regAnc = ancList[0] + "_" + ancList[1]
                    admixedIndRegions[ind][regAnc].append([regStart,regEnd,0])
                    if regEnd == Aend:
                        Apos += 1
                    else:
                        Bpos += 1
                elif Astart > Bend:
                    Bpos += 1
                elif Bstart > Aend:
                    Apos += 1
        skovInfile = "/users/kwittdil/data/data/introgression_maps/Skov_1KG/Archaicsegment_1000genomes.txt.gz"
        #skovInfile = "./Archaicsegment_1000genomes.txt.gz"
        skovSegs = {}
        for ind in popInd[pop]:
            skovSegs[ind] = []
        with gzip.open(skovInfile, 'rt') as g:
            for line in g:
                spline = line.split()
                kInd = spline[0]
                if kInd in popInd[pop]:
                    kChr = spline[1]
                    kSt = int(spline[2])
                    kEnd = int(spline[3])
                    prob = spline[7]
                    if float(prob) > 0.9 and str(kChr) == str(chromosome):
                        skovSegs[kInd].append([kSt,kEnd])
        for ancestry in diplo_ancestries:
            for ind in popInd[pop]:
                skovuRegions = skovSegs[ind]
                ancuRegions = []
                for region in admixedIndRegions[ind][ancestry]:
                    coords = region[0:2]
                    ancuRegions.append(coords)
                skovRegions = sorted(skovuRegions)
                ancRegions = sorted(ancuRegions)
                if len(skovRegions) > 0 and len(ancRegions) > 0:
                    skovRegNum = 0
                    kSt,kEnd = skovRegions[skovRegNum][0:2]
                    for regA in ancRegions:
                        aSt, aEnd = regA[0:2]
                        while int(kEnd) < int(aSt) and skovRegNum < len(skovRegions)-1: #catching up Skov to ancestry
                            skovRegNum += 1
                            kSt,kEnd = skovRegions[skovRegNum][0:2]
                        if (int(kSt) <= int(aSt) <= int(kEnd)) or (int(aSt) <= int(kSt) <= int(aEnd)): #If segments have any overlap
                            sharedLen = 0
                            mnSt = max(int(kSt),int(aSt))
                            mnEnd = min(int(kEnd), int(aEnd))
                            sharedSeg = mnEnd-mnSt
                            sharedLen += sharedSeg
                            skovRegNum += 1
                            if skovRegNum < len(skovRegions):
                                kSt,kEnd = skovRegions[skovRegNum][0:2]
                            else:
                                continue
                            while int(kSt) < int(aEnd):
                                mnSt = max(int(kSt),int(aSt))
                                mnEnd = min(int(kEnd), int(aEnd))
                                sharedSeg = mnEnd-mnSt
                                sharedLen += sharedSeg
                                skovRegNum += 1
                                if skovRegNum < len(skovRegions):
                                    kSt,kEnd = skovRegions[skovRegNum][0:2]
                                else:
                                    break
                            ancLen = int(aEnd)-int(aSt)
                            percentShared = sharedLen/ancLen
                            if percentShared < 0:
                                print(percentShared)
                            outLine = [ind, pop, ancestry, str(chromosome), str(aSt), str(aEnd), str(sharedLen), str(percentShared)]
                            tabLine = "\t".join(outLine) + "\n"
                            f.write(tabLine)   
                        elif int(aEnd) < int(kSt): #if it's not in Skov
                            regLen = int(aEnd) - int(aSt)
                            outLine = [ind, pop, ancestry, str(chromosome), str(aSt), str(aEnd), "0", "0"]
                            tabLine = "\t".join(outLine) + "\n"
                            f.write(tabLine)
                elif len(ancRegions)>0:
                    for regA in ancRegions:
                        aSt, aEnd = regA[0:2]
                        regLen = int(aEnd) - int(aSt)
                        outLine = [ind, pop, ancestry, str(chr), str(aSt), str(aEnd), "0", "0"]
                        tabLine = "\t".join(outLine) + "\n"
                        f.write(tabLine)
f.close()
