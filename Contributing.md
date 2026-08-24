# Updating EDUT (Error Detection Using Triplets) so that it can deal with modern VCF input.

## The outline of the problem
We have an older tool written in Perl called EDUT, that detected potential errors in resequencing data sets using triplets of nucleotide sites. It took both FASTA alignments and Hudson’s ms simulation format output as input. It seems like it should really now take a VCF as input. The goal would be to pass a reasonable segment (e.g., 50 SNPs) to EDUT and get output that points to individual genotypes that could be in error. One question is about phase. Genotypes could be in the 0/1 or 0|1 format. We need to determine if phasing is absolutely necessary for EDUT to make sense!
https://github.com/pmorrell/EDUT

## Additional information about EDUT
- The full text of the [Toleno et al. 2009]( https://academic.oup.com/bioinformatics/article/23/14/1807/190632?login=false) article is available from the publisher's website.
- The [EDUT manual](https://github.com/pmorrell/EDUT/blob/master/EDUT%20Manual.pdf) is in PDF format.
- The [source code](https://github.com/pmorrell/EDUT/blob/master/EDUT_1.2.1.pl) is written in Perl. Should we switch to Python3?

## EDUT can detect both issues in genotyping and in phasing
- EDUT was designed to work with phased data.
- Is the larger source of error the genotyping phase of resequencing or computational phasing?
    - For inbred lines, phasing involves a relatively small number of sites, and in many instances, it may involve only homozygous genotypes.
    - For this circumstance, where all genotypes in phased data are 0|0 or 1|1, we should consider collapsing to consider the haploid genotypes at the sites passed to edit. This would still permit tracing genotyping errors to single individuals.
    - It might be preferable to permit users to pass a VCF with no heterozygous genotypes in 0/0 or 1/1 format and run EDUT on unphased data by collapsing to haploid genotypes.
    - If a heterozygous site is encountered in unphased data, it should be reported to the user by SNP position and sample name. That genotyped position should be treated as missing data, but with the warning that "this may result in unidentified genotyping errors."

## Open questions
- What the best choice of window size? EDUT works on triplets of sites, so SNP 1 is compared to SNP 2 and 3 for the combination of three genotypes. The possible combinations are:
000 (Ancestral / All Reference)
001
010 (double recombinants or single switch errors)
011
100
101 (double recombinants or single switch errors)
110
111
- How does the algorithm scale with SNP number? It is not going to be favorable! But it may not add much information.
- In sufficiently large windows, a double recombination is possible due to shear distance between sparsely genotyped sites. We want to avoid a circumstance where errors are flagged because the SNP sites are at base positions: 1, 10000, and 100000, for example. 
- How do we handle SNP sites that overlap an indel? SNP sites could occur within the indel site. I suspect we default to counting valid genotypes, with the option to exclude them.
- How do we deal with missing data? If a site has mostly missing data, does it still go into calculations? Should users have a threshold?
- Is 50 sites too few? Perhaps compute times for 100 sites isn't a barrier? Perhaps set the window size to 100 but warn if SNPs start to become "noncontigous."
