
# ImmuneLENS

Immune Lymphocyte Estimation from Nucleotide Sequencing (Immune LENS) an R package to calculate T and B cell fractions from WGS data from hg19 or hg38 aligned genomes. Currently B cell fraction calculation is only supported for hg38 aligned genomes. It is an extension to the package T cell ExTRECT that estimated T cell content from  WES data for hg19 or hg38 aligned genomes. 

For more details on this package please read our publication [*ImmuneLENS characterizes systemic immune dysregulation in aging and cancer. Bentham et al. Nature (2025)*](https://www.nature.com/articles/s41588-025-02086-5).

## Instalation guide

Immune LENS can either be installed from github using the `install_github()` function in the package devtools or directly from source.

First make sure all dependencies are installed:
```r
dependency.packages <- c( 'dplyr', 'tidyr', 'ggnewscale', 'packcircles',
                         'readr','zoo', 'seqinr', 'rlang', 'mgcv', 'ggplot2',
                         'magrittr', 'tibble', 'devtools', 'knitr', 'restriktor',
                         'rmarkdown', 'rpart','data.table','quadprog','lavaan','ggpubr')

# Function to check and install missing packages
install_missing_packages <- function(packages) {
  installed <- installed.packages()[,"Package"]
  missing_packages <- packages[!packages %in% installed]
  
  if (length(missing_packages)) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    install.packages(missing_packages)
  } else {
    message("All packages are already installed.")
  }
}

install_missing_packages(dependency.packages)
```

### Method 1 - Using devtools

```r
# install.packages('devtools')
library(devtools)
install_github("McGranahanLab/ImmuneLENS")

```

### Method 2 - Downloaded from source

Pull ImmuneLENS repository and from terminal run `R CMD INSTALL 'PATH/TO/ImmuneLENS/'` 

Alternatively install directly from the ImmuneLENS_1.0.2.tar.gz file:

```r
install.packages('PATH/To/ImmuneLENS.tar.gz', repos=NULL, type ='source')
```

## Requirements

Samtools 

## ImmuneLENS Basic Usage 

ImmuneLENS uses a new segment based model to calculate the T cell fraction within a DNA WGS sample.             
New features include:                                                   
1. Extension to 4 V(D)J genes: TCRA, TCRB, TCRG for T cells and  IGH for B cells.                                
2. IGH class switching quantification                                
3. IGH germline and somatic haplptype calling and correction         
4. Calculation of Shannon diversity metrics based on V and J        
segment usage.                                                        
5. New visualisation functions specific to WGS samples.               
6  New function to measure distance (Jensen-Shannon) between V or J   
segment usage of two separate samples.                                  


```r
# 1. Extract coverage from BAM for TCRA, TCRB, TCRG, IGH for WGS in new function ----
library(ImmuneLENS)

# test.bam <- '/PATH/TO/BAM'
# outDir <- '/PATH/TO/OUT/DIR/'

# Note: need samtools installed and available to use
# TCRA.cov <- getCovFromBam_WGS(bamPath = test.bam, outPath = outDir,
#                               vdj.gene = 'TCRA',hg19_or_38 = 'hg38')
# 
# TCRB.cov <- getCovFromBam_WGS(bamPath = test.bam, outPath = outDir,
#                               vdj.gene = 'TCRB',hg19_or_38 = 'hg38')
# 
# TCRG.cov <- getCovFromBam_WGS(bamPath = test.bam, outPath = outDir,
#                               vdj.gene = 'TCRG',hg19_or_38 = 'hg38')
# 
# IGH.cov <- getCovFromBam_WGS(bamPath = test.bam, outPath = outDir,
#                              vdj.gene = 'IGH',hg19_or_38 = 'hg38')

# Instead for this tutorial we will use pre-extracted coverage files which are
# available on zenodo https://zenodo.org/records/11094087/

TCRA.cov <- 'https://zenodo.org/records/11094087/files/test_cov_TCRA.txt'
TCRB.cov <- 'https://zenodo.org/records/11094087/files/test_cov_TCRB.txt'
TCRG.cov <- 'https://zenodo.org/records/11094087/files/test_cov_TCRG.txt'


# Load coverage as before
TCRA.df <- loadCov(TCRA.cov)
TCRB.df <- loadCov(TCRB.cov)
TCRG.df <- loadCov(TCRG.cov)

# 2. Run T cell ExTRECT for TCRA, TCRB, TCRG  using new WGS version ----
#  runImmuneLENS
# Note - this is slower than WES version, now takes ~1-2 minutes
TCRA.out <- runImmuneLENS(vdj.region.df = TCRA.df, vdj.gene = 'TCRA',
                          hg19_or_38 = 'hg38')
TCRB.out <- runImmuneLENS(vdj.region.df = TCRB.df, vdj.gene = 'TCRB',
                          hg19_or_38 = 'hg38')
TCRG.out <- runImmuneLENS(vdj.region.df = TCRG.df, vdj.gene = 'TCRG',
                          hg19_or_38 = 'hg38')


# Output of runImmuneLENS contains a list of free objects:
#     1. A summary data frame of T/B cell fraction with calculated diversity metrics
#  of the V and J segments used according to the fitted model
TCRA.out[[1]]

#     2. A data frame containing the segment fractions for all the segments
#  fitted in the model.
TCRA.out[[2]]

#     3. A data frame containing the raw output from the model that the first two
#  objects in list were created from, contains fitted logR/T cell fraction values
TCRA.out[[3]]


# 3.If purity/ploidy is known of cancer sample can use function adjustImmuneLENS() ----
# Example:
adjustImmuneLENS(TCRA.out[[1]],purity = 0.5,
                 local.cn = 3,vdj.gene = 'TCRA')
adjustImmuneLENS(TCRA.out[[2]],purity = 0.5,
                 local.cn = 3,vdj.gene = 'TCRA')
adjustImmuneLENS(TCRA.out[[3]],purity = 0.5,
                 local.cn = 3,vdj.gene = 'TCRA')

# 3. Results can be visualised with change of segments marked by alternating colours ----
# Note this takes ~2-3 minutes
plotImmuneLENS(vdj.region.df = TCRA.df, vdj.gene = 'TCRA',
               hg19_or_38 = 'hg38', ylims = c(-1.5,1))
plotImmuneLENS(vdj.region.df = TCRB.df, vdj.gene = 'TCRB',
               hg19_or_38 = 'hg38', ylims = c(-1.5,1))
plotImmuneLENS(vdj.region.df = TCRG.df, vdj.gene = 'TCRG',
               hg19_or_38 = 'hg38', ylims = c(-1.5,1))


# 4. Alternatively to visualise diversity we can make bubble plot based on the ----
# segment proportion:
plotSegmentBubble(TCRA.out[[2]],'V',12)
plotSegmentBubble(TCRA.out[[2]],'J',12)

plotSegmentBubble(TCRB.out[[2]],'V',12)
plotSegmentBubble(TCRB.out[[2]],'J',12)

plotSegmentBubble(TCRG.out[[2]],'V',12)
plotSegmentBubble(TCRG.out[[2]],'J',12)


# 5. IGH germline haplotype calling and correction ----
IGH.germline.cov <- 'https://zenodo.org/records/11094087/files/germline_test_IGH.txt'
IGH.germline.cov.df <- loadCov(IGH.germline.cov)

germline.hap.result <- IGH_haplotype_norm_fun(IGH.germline.cov.df)
IGH.cov.df_update1 <- germline.hap.result[[1]]
IGH.regions_gl_df1 <- germline.hap.result[[2]]

# Examine the predicted non-diploid regions
IGH.regions_gl_df1

IGH.out <- runImmuneLENS(vdj.region.df = IGH.cov.df_update1, vdj.gene = 'IGH',
                                  hg19_or_38 = 'hg38',GC_correct = TRUE,
                                  removed_flag = TRUE, sample_name = 'test')

plotImmuneLENS(vdj.region.df = IGH.germline.cov.df, vdj.gene = 'IGH',
               hg19_or_38 = 'hg38', ylims = c(-1.5,1),
               sample_name = 'Without germline haplotype correction')
plotImmuneLENS(vdj.region.df = IGH.cov.df_update1, vdj.gene = 'IGH',
               hg19_or_38 = 'hg38', ylims = c(-1.5,1),
               sample_name = 'With germline haplotype correction')


# 6. IGH somatic haplotype calling and correction ---
IGH.tumour.cov <- 'https://zenodo.org/records/11094087/files/tumour_test_IGH.txt'
IGH.tumour.cov.df <- loadCov(IGH.tumour.cov)

tumour.hap.result <- IGH_haplotype_norm_fun_tumour(IGH.germline.cov.df,
                                                     IGH.tumour.cov.df,
                                                     IGH.regions_gl_df1)

IGH.tumour.cov.df_update1  <- tumour.hap.result$cov.df_update2a
IGH.tumour.cov.df_update2 <- tumour.hap.result$cov.df_update2b
somatic_regions_df <- tumour.hap.result$somatic_regions_df
somatic_correction_applied <- tumour.hap.result$somatic_correction_applied

IGH.tumour.out <- runImmuneLENS(vdj.region.df = IGH.tumour.cov.df_update2, vdj.gene = 'IGH',
                                  hg19_or_38 = 'hg38',GC_correct = TRUE,
                                  removed_flag = TRUE, sample_name = 'test')

plotImmuneLENS(vdj.region.df = IGH.tumour.cov.df, vdj.gene = 'IGH',
               hg19_or_38 = 'hg38', ylims = c(-1.5,1),
               sample_name = 'Without germline/somatic haplotype correction')
plotImmuneLENS(vdj.region.df = IGH.tumour.cov.df_update1, vdj.gene = 'IGH',
               hg19_or_38 = 'hg38', ylims = c(-1.5,1),
               sample_name = 'With germline/without somatic haplotype correction')
plotImmuneLENS(vdj.region.df = IGH.tumour.cov.df_update2, vdj.gene = 'IGH',
               hg19_or_38 = 'hg38', ylims = c(-1.5,1),
               sample_name = 'With germline + somatic haplotype correction')



plotSegmentBubble(IGH.out[[2]],'J',12)
plotSegmentBubble(IGH.tumour.out[[2]],'V',12)
```

