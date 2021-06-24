# Prerequisites

R packages: 
- PMA
- mvtnorm
- GBJ
- TisCoMM
- Rcpp
- plink2R (https://github.com/gabraham/plink2R)

Other packages: for preprocess of genotype data
- [Tabix](http://www.htslib.org/doc/tabix.html) 
- [plink2](https://www.cog-genomics.org/plink/2.0/) 

# Simulations

In the scripts/simulations/ directory, there are scripts we used to conduct 1) Type I error evaluation, and 2) Power analysis for our MSG approach and S-MultiXcan, UTMOST, and sCCA+ACAT. 

# Real data analysis

In the scripts/real_data/ directory, there are scripts we used to apply our MSG approach and established S-MultiXcan, UTMOST, and sCCA+ACAT to real data from GTEx

The application of MSG to real data consists of the following steps

## Generate files for X and Y matrix

scripts/real_data/step0_format_gtex_data_cmds.pl and scripts/real_data/step1_genMatrix_gtex.pl are Perl scripts written by Qiang Wei to process the genotype and splicing event data provided on GTEx Protal. The output of these scripts are genotype matrix (X) and splicing event matrix (Y matrix). Examples can be found in example_data/ folder.

## Train MSG models

We applied sCCA to the X,Y matrices generated in the previous step

## Summarize and plots

In this step, we summarize the results from all tested genes for a trait, and produce 1) bar plot, 2) venn diagram, 3) manhatten plots to visualize the results by different methods.

# Acknowledgement

Part of the code is modified from previous work: 
- S-MultiXcan at https://github.com/hakyimlab/MetaXcan
- TWAS at https://github.com/gusevlab/fusion_twas/
- UTMOST at https://github.com/Joker-Jerome/UTMOST
- TisCoMM at https://github.com/XingjieShi/TisCoMM
- JTI at https://github.com/gamazonlab/MR-JTI

We thank the authors for sharing their code! 

