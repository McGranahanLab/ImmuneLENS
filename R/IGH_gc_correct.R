#' GC correct IGH coverage values for use prior to identification of IGH haplotype
#'
#' @param input_cov data frame of positions and coverage values for IGH loci
#' @param hg19_or38 default = 'hg38'
#' @return  Data frame of GC corrected coverage values
#' @name IGH_gc_correct

IGH_gc_correct <- function(input_cov, hg19_or38 = 'hg38'){
  X2 <- pos <- GC <- exon.GC <- smooth.gc <- reads.gc.correct <- reads <- exon.gc <-  NULL
  
  # Run GC correction on coverage values before running haplotype norm
  # Use constrained linear model (restriktor)
  if(hg19_or38 == 'hg38'){
    VDJ_fasta <- ImmuneLENS::all_fasta_hg38[['IGH']]
    vdj.seg <- vdj_seg_list[['IGH_hg38']]
  }else{
    VDJ_fasta <- ImmuneLENS::all_fasta[['IGH']]
    vdj.seg <- vdj_seg_list[['IGH_hg19']]
  }
  vdj.chr <- 'chr14'
  vdj.start <- vdj.seg[1,2]
  vdj.end <- vdj.seg[1,3]
  
  IGH.exons <- data.frame(X1 = vdj.chr,
                          X2 = seq(from = vdj.start, to = vdj.end - 1000,
                                   by = 1000)) %>%
    dplyr::mutate(X3 = X2 + 999)
  vdj.end.but1 <- IGH.exons$X3[dim(IGH.exons)[1]] + 1
  IGH.exons <- rbind(IGH.exons,
                     data.frame(X1 = vdj.chr, X2 = vdj.end.but1, X3 = vdj.end))
  
  # Get GC content for exons
  exon.adjust.loc <- vdj.start -1
  exons.selected <- IGH.exons
  IGH.exons.loc <- list()
  for(i in seq_len(dim(exons.selected)[1])){
    IGH.exons.loc[[i]] <- c(exons.selected$X2[i] - exon.adjust.loc,
                            exons.selected$X3[i] - exon.adjust.loc)
  }
  exons.gc.content <- exonwindowplot2(IGH.exons.loc, VDJ_fasta[[1]],0)
  
  gc.df <- slidingwindowplot_alt(1000, VDJ_fasta[[1]])
  gc.df$pos <- gc.df$loc + exon.adjust.loc
  gam.model <- mgcv::gam(GC~s(pos, bs = 'cs'), data = gc.df)
  
  get_gc_prediction <- function(x){
    mgcv::predict.gam(gam.model, newdata = data.frame(pos = x))}
  
  input_cov_gc <- input_cov %>%
    dplyr::mutate(exon = exonPosFun_v(pos, IGH.exons)) %>%
    dplyr::left_join(exons.gc.content, 'exon') %>%
    dplyr::rename(exon.gc = GC) %>%
    dplyr::mutate(exon.gc2 = exon.gc^2) %>%
    dplyr::mutate(smooth.gc = get_gc_prediction(pos)) %>%
    dplyr::mutate(smooth.gc2 = smooth.gc^2)
  
  # Finally normalise by the GC content
  gc.lm <- lm(reads ~ exon.gc + exon.gc2 + smooth.gc + smooth.gc2,
              y = TRUE, data = input_cov_gc)
  # Do constrained linear model using restriktor - what are the constraints???? 
  # Running test on data to see:
  myConstraints <- 'smooth.gc > 10 \n smooth.gc2 < 0 \n -1*smooth.gc2 == smooth.gc'
  restr.lm <- restriktor::conLM.lm(gc.lm,
                                   constraints = myConstraints,
                                   se = 'none',mix.weights = 'none',
                                   control = list(absval = 1e-5))
  if(is.null(dim(restr.lm$residuals))){
    res.df <- restr.lm$residual
  }else{
    res.df <- restr.lm$residuals[,1]
  }
  
  gc.lm.df <- data.frame(pos = input_cov_gc$pos[as.numeric(names(res.df))],
                         reads.gc.correct = res.df)  
  
  
  
  if(min(gc.lm.df$reads.gc.correct) < 0){
    gc.lm.df$reads.gc.correct <- gc.lm.df$reads.gc.correct + abs(min(gc.lm.df$reads.gc.correct))
  }
  
  cov.df.gc.update <- gc.lm.df %>%
    dplyr::rename(reads = reads.gc.correct)
  # Finally make sure old zeros are still zeros
  zero_pos <- input_cov %>% 
    dplyr::filter(reads == 0) %>% 
    dplyr::select(pos) %>% `[[`(1)
  zero_peak <- cov.df.gc.update %>% 
    dplyr::filter(pos %in% zero_pos) %>% dplyr::select(reads) %>% `[[`(1) %>% 
    median(na.rm = TRUE)
  
  cov.df.gc.update$reads <- cov.df.gc.update$reads - zero_peak
  cov.df.gc.update <- cov.df.gc.update %>% 
    dplyr::mutate(reads = ifelse(reads < 0, 0, reads))
  
  
  return(cov.df.gc.update)
}