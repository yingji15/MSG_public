This directory includes our example code and data for manuscript "Integration of multidimensional splicing data and GWAS summary statistics for risk gene discovery"

# scripts/

This folder includes example scripts to conduct the analysis in our manuscript.

## Prerequisites

Please install the following packages before running the scripts, as the scripts depend on some functions from these packages

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


## Simulations

In the scripts/simulations/ directory, there are scripts we used to conduct 1) Type I error evaluation, and 2) Power analysis for our MSG approach and S-MultiXcan, UTMOST, and sCCA+ACAT. 

## Real data analysis

In the scripts/real_data/ directory, there are scripts we used to apply our MSG approach and established S-MultiXcan, UTMOST, and sCCA+ACAT to real data from GTEx

The application of MSG to real data consists of the following steps

### Generate files for X and Y matrix

scripts/real_data/step0_format_gtex_data_cmds.pl and scripts/real_data/step1_genMatrix_gtex.pl are Perl scripts written by Qiang Wei to process the genotype and splicing event data provided on GTEx Protal. 

Input: 

- genotype: vcf.gz file
- splicing phenotype: gtex processed splicing data from leafcutter, in bed.gz format, example: https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_sQTL_phenotype_matrices.tar
- files with reference data: provided in ref_data/, here we provided hg38 version 
- 
Output: genotype matrix (X) and splicing event matrix (Y matrix). Examples can be found in example_data/ folder.

### Train MSG models

In scripts/real_data/2_MSG_models.R, we applied sCCA to the X,Y matrices generated in the previous step

Input: X and Y matrices generated in previous step

Output: RData file with weights


### Summarize and plots

In scripts/real_data/step3_summarize_and_plots.R, we summarize the results from all tested genes for a trait, and produce 1) bar plot, 2) venn diagram, 3) manhatten plots to visualize the results by different methods.


# example_data/

As mentioned earlier, this folder includes several example datasets for real data application

# ref_data/

This folder includes some public available datasets on coding genes, with information of their location etc. These information are used in the preprocess of genotype data, and generate LD panels.  

# Acknowledgement

Part of the code is modified from previous work: 
- S-MultiXcan at https://github.com/hakyimlab/MetaXcan
- TWAS at https://github.com/gusevlab/fusion_twas/
- UTMOST at https://github.com/Joker-Jerome/UTMOST
- TisCoMM at https://github.com/XingjieShi/TisCoMM
- JTI at https://github.com/gamazonlab/MR-JTI

We thank the authors for sharing their code! 

