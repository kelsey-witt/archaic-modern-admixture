# archaic-modern-admixture
This repository contains the code needed to repeat the analyses used in "The impact of modern admixture on archaic human ancestry
in human populations", published in GBE in 2023 by Witt and colleagues. A link to the publication will be added once the manuscript is officially published.

The aim of this study was to examine how archaic ancestry is distributed across modern ancestry segments in admixed populations. For this study we specifically used the admixed American populations from the 1000 Genomes dataset, archaic ancestry calls using Sprime from Browning et al. 2018 (https://www.cell.com/fulltext/S0092-8674(18)30175-2), modern ancestry tracts from Martin et al. 2017 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5384097/), and ancestry tract data from Sankararaman et al. 2016 (https://www.nature.com/articles/nature12961) and Skov et al. 2018 (https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007641).

## Data Used
* [1000 Genomes Project Phase 3 VCF files](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/)
* [Modern Ancestry Tract Calls for Admixed American 1000 Genomes Populations](https://personal.broadinstitute.org/armartin/tgp_admixture/)
* [Sprime archaic site calls](https://data.mendeley.com/datasets/y7hyt83vxr/1)

## Analyses and scripts
### Calculating archaic allele density
archaic_snp_density_sprime_admixed.py divides the genomes of the admixed American populations (currently set to PEL, MXL, CLM, and PUR) into ancestry segments, and then calculates the % of sites within those ancestry segments that are classified as archaic by Sprime.
* Usage: python3 archaic_snp_density_sprime_admixed.py [archaic_set], where archaic_set refers to the set of archaic sites under analysis. Options for archaic_set are nd_either, nd_both, n_only, n_all, d_only, d_all, and are explained fully at the top of the script.
* Inputs: 1000 Genomes vcf and panel file, Sprime calls, and modern ancestry tract calls. Input file paths are all specified in the script itself and will need to be altered to fit your organization system. (Variables that will need to be checked are pop_file (1000 Genomes panel file, line 26), outfile (output file name, line 25), infile (the ancestry tract bed files, line 83), modern_file (vcf file, line 137), and sprime_file (sprime calls, line 140)
* Output: A csv file where each row is a different ancestry tract in a different individual and the columns are, in order, the population and ancestry call separated by an underscore, the individual ID, the coordinates of the tract in standard BED format, the number of archaic alleles, the length of the tract, and the density (the number of archaic alleles divided by the tract length)

### Counting archaic SNPs per individual
archaic_snps_perind_sprime.py counts the number of sites with archaic SNPs, and the total number of archaic SNPs (ie a homozygous site would count as 2 SNPs) for each individual, using the Sprime calls to define archaic sites.
* Usage: python3 archaic_snps_perind_sprime.py [archaic_set], where archaic_set refers to the set of archaic sites under analysis. Options for archaic_set are nd_either, nd_both, n_only, n_all, d_only, d_all, and are explained fully at the top of the script.
* Inputs: 1000 Genomes vcf and panel file, and Sprime calls. Input file paths are all specified in the script itself and will need to be altered to fit your organization system. (Variables that will need to be checked are pop_file (1000 Genomes panel file, line 24), outfile (output file name, line 25), modern_file (vcf file, line 60), and sprime_file (sprime calls, line 64)
* Output: A csv file where each row is a different individual and the columns are, in order, the individual, their population, the number of sites with archaic SNPs and the total number of archaic SNPs across the genome

### Calculating heterozygosity for ancestry tracts
calc_heterozygosity_per_ancestry_v2.py calculates the heterozygosity for each ancestry tract in each individual in the admixed American populations (PEL, PUR, CLM, MXL).
* Usage: python3 calc_heterozygosity_per_ancestry_v2.py
* Inputs: 1000 Genomes vcf and panel file, and modern ancestry tract calls. Input file paths are all specified in the script itself and will need to be altered to fit your organization system. (Variables that will need to be checked are pop_file (1000 Genomes panel file, line 19), outfile (output file name, line 18), infile (the ancestry tract bed files, line 90), and modern_file (vcf file, line 143)
* Outputs: A csv file where each row is a different individual and ancestry combination and the columns are, in order, population, ancestry, individual ID, chromosome, region number (an internal counter of ancestry tracts), position start position end, and number of heterozygous sites

