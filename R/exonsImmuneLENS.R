#' Function to run ImmuneLENS for exon output
#'
#' @param vdj.region.df data frame containing coverage values by position
#' @param vdj.gene VDJ gene to use
#' @param hg19_or_38 hg19 or hg38 version of genome
#' @param GC_correct whether to use GC correction or not for output
#' @param median.k rolling median window
#' @param median.thresh threshold to remove exons with low coverage
#' @param sample_name name of sample run
#' @param customSeg custom seg file
#' @param customFasta custom Fasta file
#' @param customFlag custom list of exons to flag for removal
#' @param removed_flag Whether to removed flagged WGS locations (default = TRUE)
#' @return data frame of TCRA T cell fractions with 95\% CI
#' @name exonsImmuneLENS
#' @export
exonsImmuneLENS <- function(vdj.region.df, vdj.gene = 'TCRA',
                                         hg19_or_38 = 'hg38',
                                         GC_correct = TRUE,
                                         median.k = 50, median.thresh = 15,
                                         sample_name = 'test',
                                         customSeg = NULL,
                                         customFasta = NULL,
                                         customFlag = NULL,
                                         removed_flag = TRUE){

  X2 <- exon2 <- pos <-  NULL
  # For WGS no exons are required
  # colnames(exons.selected) <- c('X1','X2','X3')
  vdj.chr.df <- data.frame(gene = c('TCRA','TCRB','TCRG','IGH','IGL','IGK', 'TCRD'),
                           chr = c('chr14','chr7','chr7','chr14','chr22','chr2','chr14'))

  if(is.null(customSeg)){
    # data("vdj_seg_list")
    seg.name <- paste0(vdj.gene, '_', hg19_or_38)
    vdj.seg <- vdj_seg_list[[seg.name]]
  }else{
    vdj.seg <- customSeg
  }


  if(!is.null(customFasta)){
    VDJ_fasta <- seqinr::read.fasta(customFasta)
  }else{
    if(hg19_or_38 == 'hg19'){
      # data("all_fasta")
      VDJ_fasta <- all_fasta[[vdj.gene]]
    }else{
      # data("all_fasta_hg38")
      VDJ_fasta <- all_fasta_hg38[[vdj.gene]]
    }
  }


  vdj.chr <- vdj.chr.df$chr[which(vdj.chr.df$gene == vdj.gene)]
  vdj.start <- vdj.seg[1,2]
  vdj.end <- vdj.seg[1,3]

  TCRA.exons <- data.frame(X1 = vdj.chr,
                           X2 = seq(from = vdj.start,
                                    to = vdj.end - 1000, by = 1000)) %>%
    dplyr::mutate(X3 = X2 + 999)
  vdj.end.but1 <- TCRA.exons$X3[dim(TCRA.exons)[1]] + 1

  TCRA.exons <- rbind(TCRA.exons,
                      data.frame(X1 = vdj.chr, X2 = vdj.end.but1, X3 = vdj.end))

  exons.selected <- TCRA.exons

  # Get GC content for exons
  exon.adjust.loc <- vdj.start - 1

  TCRA.exons.loc <- list()
  for(i in seq_len(dim(exons.selected)[1])){
    TCRA.exons.loc[[i]] <- c(exons.selected$X2[i] - exon.adjust.loc,
                             exons.selected$X3[i] - exon.adjust.loc)

  }
  # These are the GC content within the exons


  exons.gc.content <- exonwindowplot2(TCRA.exons.loc, VDJ_fasta[[1]],0)

  if(removed_flag){
    if(!is.null(customFlag)){
      exons.flaged <- customFlag
    }else{
      exons.flaged <- flagged_exons[[vdj.gene]][[1]]
    }
    vdj.region.df <-vdj.region.df %>%
      dplyr::mutate(exon2 = (pos - vdj.start + 100) %/% 100) %>%
      dplyr::filter(!exon2 %in% exons.flaged) %>%
      dplyr::select(-exon2)
  }

  # 1. Check VDJ file:
  if(dim(vdj.region.df)[1] == 0){
    return(NULL)
  }

  exons.to.use <- seq_len(dim(TCRA.exons)[[1]])
  # median filter coverage for normalisation
  median.exon.output <- medianExonCoverage(vdj.region.df, exons.selected,
                                           median.k, median.thresh,
                                           exons.to.use)
  vdj.region.df.filt.exons.median <- median.exon.output[[1]]
  exon.remove <- median.exon.output[[2]]
  
  if(dim(vdj.region.df.filt.exons.median)[1] == 0) stop('All positions have been removed due to low coverage. Consider decreasing median.thresh')
  

  # Calculate log ratio
  vdj.logR.df <- getLogRdf(vdj.region.df.filt.exons.median, vdj.seg, minCov = 0)

  # VDJ.QC check
  vdj.logR.df <- vdj.logR.df[!is.infinite(vdj.logR.df$Ratio), ]
  vdj.logR.df <- vdj.logR.df[!is.na(vdj.logR.df$Ratio), ]

  if(dim(vdj.logR.df)[1] == 0 | length(exon.remove) > 30000){
    return(NA)}

  if(GC_correct){
    vdj.logR.df <- GCcorrect_WGS(vdj.logR.df,
                                 exons = exons.selected,
                                 exonList = exons.gc.content,
                                 gene.fasta = VDJ_fasta,
                                 gene.fasta.start = exon.adjust.loc,
                                 hg19_or38 = hg19_or_38)
    baselineAdj.out <- baselineAdj(vdj.logR.df,
                                   vdj.seg, GCcorrect = TRUE)
    vdj.logR.df <- baselineAdj.out[[1]] 

   # vdj.fraction.output  <-  getVDJfraction_segmodel(vdj.logR.df, sample_name, hg19_or_38)

  }else{
    baselineAdj.out <- baselineAdj(vdj.logR.df, vdj.seg, GCcorrect = FALSE)
    vdj.logR.df <-baselineAdj.out[[1]]

   # vdj.fraction.output  <- getVDJfraction_segmodel(vdj.logR.df, sample_name,hg19_or_38)
  }

  #out.df <- vdj.fraction.output
  return(vdj.logR.df)
}
