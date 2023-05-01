"""
This script was written by Kelsey Witt Dillon in February 2023. The goal of this script is to 
identify how archaic introgression tracts from Sankararaman et al. 2016 overlap with modern
ancestry tracts in admixed American individuals from Martin et al. 2017. The haploid ancestry tracts
from Sankararaman are combined into diploid ancestry tracts - if one chromosome has archaic ancestry,
the region is considered archaic. The genomes of the individuals are subdivided into regions of 
diploid ancestry (ie AFR_AFR, EUR_AFR, etc.) and then the percentage of that ancestry tract that 
overlaps with a diploid archaic ancestry tract is calculated.

The input files are the .panel file from the 1000 genomes dataset ("pop_file", line 25),
the ancestry bed files from Martin et al., ("infile", line 79), the population ID file
from Sankararaman ("sankPopFile", line 117), and the archaic haplotype data ("sankInfile", line 126)

Usage: python3 sank_region_overlap_AAA_duoAnc.py

This script generates a tab-delimited outfile, currently named "sank_ancestry_overlap.txt" ("outfile", 
line 66). The columns are the individual, the population, the diploid ancestry of the region, the 
chromosome, the start of the ancestry region, the end of the ancestry region, the length of the region
that includes archaic ancestry, and the percentage that length represents.
"""

import gzip
import sys

pop_file = "/users/kwittdil/data/data/modern_genomes/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
admixPops = ["PUR", "CLM", "MXL"]
unrepSankInd = ["HG01278", "HG01274", "NA19737", "NA19660", "NA19685", "NA19753", "NA19675", "NA19738", "NA19672"] #some sank individuals aren't in panel file

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

outputFile = "./sank_ancestry_overlap.txt"
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

        sankSegs = {}
        for ind in popInd[pop]:
            sankSegs[ind] = []
        sankPopFile = "/users/kwittdil/data/data/introgression_maps/Sankararaman_1KG_ALL/package/ids/" + pop + ".ids"
        sankIndOrder = []
        with open(sankPopFile) as d:
            indName = "null"
            for line in d:
                dspline = line.split()
                dID = dspline[0]
                if dID not in sankIndOrder: # only adds them once
                    sankIndOrder.append(dID)
        sankInfile = "/users/kwittdil/data/data/introgression_maps/Sankararaman_1KG_ALL/package/" + pop + ".hapmap/summaries/haplotypes/chr-" + chromosome + ".thresh-90.length-0.00.haplotypes"
        with open(sankInfile) as h:
            for line in h:
                if "##" not in line:
                    hspline = line.split()
                    indNum, aSt, aEnd = hspline[1:4]
                    indPos = int(indNum)//2 #divided by 2 no remainder
                    indID = sankIndOrder[indPos]
                    if indID not in unrepSankInd:
                        sankSegs[indID].append([aSt,aEnd])
        #sorting sank file to be in order
        sankOrder = []
        for sankid in sankIndOrder:
            if sankid not in unrepSankInd:
                sankOrder.append(sankid)
        for ind in sankOrder:
            sortSegments = sorted(sankSegs[ind])
            summedRegions = []
            if len(sortSegments)>0:
                regionNum = 0
                currSt, currEnd = sortSegments[regionNum][0:2]
                while regionNum < len(sortSegments)-1:
                    nextSt, nextEnd = sortSegments[regionNum+1][0:2]
                    if int(currEnd) < int(nextSt):
                        summedRegions.append([str(currSt),str(currEnd)])
                        regionNum += 1
                        currSt, currEnd = sortSegments[regionNum][0:2]
                    elif int(currSt) == int(nextSt): #if one region is completely contained in the next, move on
                        regionNum += 1
                        currSt, currEnd = sortSegments[regionNum][0:2]
                    elif int(nextEnd) <= int(currEnd): #if next region is contained in the previous, retain current and move on
                        regionNum += 1
                    elif int(currEnd)>=int(nextSt):
                        currEnd = max(currEnd,nextEnd)
                        regionNum += 1
                summedRegions.append([str(currSt),str(currEnd)])
            sankSegs[ind] = summedRegions
        for ancestry in diplo_ancestries:
            for ind in popInd[pop]:
                sankuRegions = sankSegs[ind]
                ancuRegions = []
                for region in admixedIndRegions[ind][ancestry]:
                    coords = region[0:2]
                    ancuRegions.append(coords)
                sankRegions = sorted(sankuRegions)
                ancRegions = sorted(ancuRegions)
                if len(sankRegions) > 0 and len(ancRegions) > 0:
                    sankRegNum = 0
                    kSt,kEnd = sankRegions[sankRegNum][0:2]
                    for regA in ancRegions:
                        aSt, aEnd = regA[0:2]
                        while int(kEnd) < int(aSt) and sankRegNum < len(sankRegions)-1: #catching up Sank to ancestry
                            sankRegNum += 1
                            kSt,kEnd = sankRegions[sankRegNum][0:2]
                        if (int(kSt) <= int(aSt) <= int(kEnd)) or (int(aSt) <= int(kSt) <= int(aEnd)): #If segments have any overlap
                            sharedLen = 0
                            mnSt = max(int(kSt),int(aSt))
                            mnEnd = min(int(kEnd), int(aEnd))
                            sharedSeg = mnEnd-mnSt
                            sharedLen += sharedSeg
                            sankRegNum += 1
                            if sankRegNum < len(sankRegions):
                                kSt,kEnd = sankRegions[sankRegNum][0:2]
                            else:
                                continue
                            while int(kSt) < int(aEnd):
                                mnSt = max(int(kSt),int(aSt))
                                mnEnd = min(int(kEnd), int(aEnd))
                                sharedSeg = mnEnd-mnSt
                                sharedLen += sharedSeg
                                sankRegNum += 1
                                if sankRegNum < len(sankRegions):
                                    kSt,kEnd = sankRegions[sankRegNum][0:2]
                                else:
                                    break
                            ancLen = int(aEnd)-int(aSt)
                            percentShared = sharedLen/ancLen
                            if percentShared < 0:
                                print(percentShared)
                            outLine = [ind, pop, ancestry, str(chromosome), str(aSt), str(aEnd), str(sharedLen), str(percentShared)]
                            tabLine = "\t".join(outLine) + "\n"
                            f.write(tabLine)   
                        elif int(aEnd) < int(kSt): #if it's not in Sank
                            regLen = int(aEnd) - int(aSt)
                            outLine = [ind, pop, ancestry, str(chromosome), str(aSt), str(aEnd), "0", "0"]
                            tabLine = "\t".join(outLine) + "\n"
                            f.write(tabLine)
                elif len(ancRegions)>0:
                    for regA in ancRegions:
                        aSt, aEnd = regA[0:2]
                        regLen = int(aEnd) - int(aSt)
                        outLine = [ind, pop, ancestry, str(chromosome), str(aSt), str(aEnd), "0", "0"]
                        tabLine = "\t".join(outLine) + "\n"
                        f.write(tabLine)
f.close()
