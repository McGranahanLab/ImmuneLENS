
# ImmuneLENS

Immune Lymphocyte Estimation from Nucleotide Sequencing (Immune LENS) an R package to calculate T and B cell fractions from WGS data from hg19 or hg38 aligned genomes. It also includes functions related to the measure of T cell fractions from WES data for hg19 or hg38 aligned genomes originally used in the R package T cell ExTRECT. For more details on the WES version please read our publication [*Using DNA sequencing data to quantify T cell fraction and therapy response. Bentham et al. Nature (2021)*](https://www.nature.com/articles/s41586-021-03894-5).

## Instalation guide

Immune LENS can either be installed from github using the `install_github()` function in the package devtools or directly from source.

### Method 1 - Using devtools

```r
# install.packages('devtools')
library(devtools)
install_github("McGranahanLab/ImmuneLENS")

```

### Method 2 - Downloaded from source


```r
install.packages('PATH/To/ImmuneLENS/', repos=NULL, type ='source')
```

## Requirements

Samtools (>v1.3.1)

## Example use
Running ImmuneLENS on your data is both fast and easy!

```r
library(ImmuneLENS)
```

First take a WGS  bam file that has been aligned to either hg19 or hg38
```r
test.bam <- '/PATH/TO/BAM/'
outDir <- '/PATH/TO/OUTPUT/DIR/'
```
Generate coverage data for the TCRA, TCRB, TCRG and IGH loci using the `getCovFromBam_WGS` function

```r
TCRA.cov <- getCovFromBam_WGS(bamPath = test.bam, outPath = outDir,
                              vdj.gene = 'TCRA',hg19_or_38 = 'hg38')
TCRB.cov <- getCovFromBam_WGS(bamPath = test.bam, outPath = outDir,
                              vdj.gene = 'TCRB',hg19_or_38 = 'hg38')
TCRG.cov <- getCovFromBam_WGS(bamPath = test.bam, outPath = outDir,
                              vdj.gene = 'TCRG',hg19_or_38 = 'hg38')
IGH.cov <- getCovFromBam_WGS(bamPath = test.bam, outPath = outDir,
                              vdj.gene = 'IGH',hg19_or_38 = 'hg38')
```

Load coverage as before 

```r
TCRA.df <- loadCov(TCRA.cov)
TCRB.df <- loadCov(TCRB.cov)
TCRG.df <- loadCov(TCRG.cov)
IGH.df <- loadCov(IGH.cov)
```

Run ImmuneLENS (note this takes 2-3 minutes per sample)

```r
TCRA.out <- runImmuneLENS(vdj.region.df = TCRA.df, vdj.gene = 'TCRA',
                                         hg19_or_38 = 'hg38',GC_correct = TRUE,
                                         removed_flag = TRUE)
TCRB.out <- runImmuneLENS(vdj.region.df = TCRB.df, vdj.gene = 'TCRB',
                                         hg19_or_38 = 'hg38',GC_correct = TRUE,
                                         removed_flag = TRUE)
TCRG.out <- runImmuneLENS(vdj.region.df = TCRG.df, vdj.gene = 'TCRG',
                                         hg19_or_38 = 'hg38',GC_correct = TRUE,
                                         removed_flag = TRUE)
IGH.out <- runImmuneLENS(vdj.region.df = IGH.df, vdj.gene = 'IGH',
                                         hg19_or_38 = 'hg38',GC_correct = TRUE,
                                        removed_flag = TRUE)
```

Output of runTcellExTRECT_WGS_segmodel contains a list of three objects:

1. A summary data frame of T/B cell fraction with calculated diversity metrics of the V and J segments used according to the fitted model
 ```r
TCRA.out[[1]]
```
2. A data frame containing the segment fractions for all the segments fitted in the model.
 ```r
TCRA.out[[2]]
```
3. A data frame containing the raw output from the model that the first two  objects in list were created from, contains fitted logR/T cell fraction values and model R2, log liklihood and S2 values
 ```r
TCRA.out[[3]]
```


If purity/ploidy is known of cancer sample can use function `adjustImmuneLENS()`

```r
adjustTcellExTRECT_WGS(TCRA.out[[1]],purity = 0.5,
                       local.cn = 3,vdj.gene = 'TCRA')
adjustTcellExTRECT_WGS(TCRA.out[[2]],purity = 0.5,
                       local.cn = 3,vdj.gene = 'TCRA')
adjustTcellExTRECT_WGS(TCRA.out[[3]],purity = 0.5,
                       local.cn = 3,vdj.gene = 'TCRA')
```

Results can be visualised with change of segments marked by alternating colours.
Note this is also slow as is running the model independently of runImmuneLENS

```r
plotImmuneLENS(vdj.region.df = TCRA.df, vdj.gene = 'TCRA',
                                 hg19_or_38 = 'hg38', ylims = c(-1.5,1))
plotImmuneLENS(vdj.region.df = TCRB.df, vdj.gene = 'TCRB',
                     hg19_or_38 = 'hg38', ylims = c(-1.5,1))
plotImmuneLENS(vdj.region.df = TCRG.df, vdj.gene = 'TCRG',
                     hg19_or_38 = 'hg38', ylims = c(-1.5,1))
plotImmuneLENS(vdj.region.df = IGH.df, vdj.gene = 'IGH',
                     hg19_or_38 = 'hg38', ylims = c(-1.5,1))
```
Alternatively to visualise diversity we can make bubble plot based on the segment proportion:

```r
plotSegmentBubble(TCRA.out[[2]],'V',12)
plotSegmentBubble(TCRA.out[[2]],'J',12)

plotSegmentBubble(TCRB.out[[2]],'V',12)
plotSegmentBubble(TCRB.out[[2]],'J',12)

plotSegmentBubble(TCRG.out[[2]],'V',12)
plotSegmentBubble(TCRG.out[[2]],'J',12)

plotSegmentBubble(IGH.out[[2]],'V',12)
plotSegmentBubble(IGH.out[[2]],'J',12)
```

To compare the segment output of two solutions can use the JSD function

```r
test.bam2 <- '/PATH/TO/EXAMPLE/BAM2'
TCRA.cov2 <- getCovFromBam_WGS(bamPath = test.bam2, outPath = outDir,
                              vdj.gene = 'TCRA',hg19_or_38 = 'hg38')

TCRA.df2 <- loadCov(TCRA.cov2)
TCRA.out2 <- runTcellExTRECT_WGS_segmodel(vdj.region.df = TCRA.df2, vdj.gene = 'TCRA',
                                         hg19_or_38 = 'hg38',GC_correct = TRUE,
                                         removed_flag = TRUE)


calcJSD_distance(TCRA.out[[2]], TCRA.out2[[2]],'V')
calcJSD_distance(TCRA.out[[2]], TCRA.out2[[2]],'J')
```
