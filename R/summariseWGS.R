#' Function to summarise raw output from WGS version
#'
#' @param segmodel.out output from getVDJfraction_WGS_segmodel
#' @param vdj.gene VDJ gene to use
#' @param customNorm Custom dataframe with average values to normalise for, summarised over 100bp (default = NULL)
#' @importFrom rlang :=
#' @name summariseWGS

summariseWGS <- function(segmodel.out, vdj.gene = 'TCRA', customNorm = NULL){
  segment <- logR <- modelR2 <- model.loglik <- model.s2 <- segment.fraction <- NULL
  model.BIC <- NULL
  p <- shannon.value <- NULL

  if(vdj.gene %in% c('TCRA','TCRB','TCRG')){
    fraction_col <- rlang::sym(paste0(vdj.gene,'.tcell.fraction'))
    vseg <- paste0(gsub('TCR','TR',vdj.gene),'V')
    jseg <- paste0(gsub('TCR','TR',vdj.gene),'J')
  }else if(vdj.gene == 'IGH'){
    fraction_col <- rlang::sym(paste0(vdj.gene,'.bcell.fraction'))
    vseg <- 'IGHV'
    jseg <- 'IGHJ'
  }
  focal_segment <- segmodel.out$segment[grepl(vseg, segmodel.out$segment) &
                                          grepl(jseg, segmodel.out$segment)]
  if(length(focal_segment) == 0){
    warning('No focal segment found - assumed low quality sample - use last vseg if possible')
    focal_segment <- segmodel.out$segment[grepl(vseg, segmodel.out$segment)]
    if(length(focal_segment) == 0){
      return(NULL)
    }else{
      focal_segment <- rev(focal_segment)[1]
      warning(paste0('focal segment set to ', focal_segment))
    }
  }
  focal_segment_split <- strsplit(focal_segment,'_')[[1]]

  output.summary <- segmodel.out %>%
    dplyr::filter(segment %in% focal_segment) %>%
    dplyr::select(sample,logR, !!fraction_col, modelR2, model.loglik, model.s2, model.BIC)

  output.segments <- segmodel.out %>%
    dplyr::filter(!segment %in% c('exon.gc','exon.gc2','smooth.gc','smooth.gc2'))

  if(!is.null(customNorm)){
    norm.custom.cols <- setdiff(colnames(customNorm),'exon2')
    segmodel.out %>%
      dplyr::filter(!segment %in% norm.custom.cols)
  }


  # 1. Get summary values of model
  if(vdj.gene %in% c('TCRA','TCRB')){
    output.v.segments <- output.segments %>%
      dplyr::filter(grepl(vseg,segment)) %>%
      dplyr::mutate(segment.fraction = !!fraction_col - dplyr::lag(!!fraction_col, default = 0)) %>%
      dplyr::mutate(segment.fraction = ifelse(segment.fraction < 0, 0, segment.fraction)) %>%
      dplyr::mutate(segment = gsub(paste0('_',focal_segment_split[2]),'',segment)) %>%
      dplyr::select(sample, segment, segment.fraction)

    output.j.segments <- output.segments %>%
      dplyr::filter(grepl(jseg, segment)) %>%
      dplyr::mutate(segment.fraction = !!fraction_col - dplyr::lead(!!fraction_col,default = 0)) %>%
      dplyr::mutate(segment.fraction = ifelse(segment.fraction < 0, 0, segment.fraction)) %>%
      dplyr::mutate(segment = gsub(paste0(focal_segment_split[1],'_'),'',segment)) %>%
      dplyr::select(sample, segment, segment.fraction)

    if(vdj.gene == 'TCRB'){
      # Fix TCRB manually for now
      output.v.segments$segment[length(output.v.segments$segment)] <- 'TRBV29_1'
      output.j.segments$segment[1] <- 'TRBJ1_1'
    }

    output.all.segments <- rbind(output.v.segments, output.j.segments)

  }else if(vdj.gene %in% c('TCRG','IGH')){
    output.j.segments <- output.segments %>%
      dplyr::filter(grepl(jseg,segment)) %>%
      dplyr::mutate(segment.fraction = !!fraction_col - dplyr::lag(!!fraction_col, default = 0)) %>%
      dplyr::mutate(segment.fraction = ifelse(segment.fraction < 0, 0, segment.fraction)) %>%
      dplyr::mutate(segment = gsub(paste0('_',focal_segment_split[2]),'',segment)) %>%
      dplyr::select(sample, segment, segment.fraction)

    output.v.segments <- output.segments %>%
      dplyr::filter(grepl(vseg, segment)) %>%
      dplyr::mutate(segment.fraction = !!fraction_col - dplyr::lead(!!fraction_col,default = 0)) %>%
      dplyr::mutate(segment.fraction = ifelse(segment.fraction < 0, 0, segment.fraction)) %>%
      dplyr::mutate(segment = gsub(paste0(focal_segment_split[1],'_'),'',segment)) %>%
      dplyr::select(sample, segment, segment.fraction)

    output.all.segments <- rbind(output.j.segments, output.v.segments)
  }

  # Add in class switching for IGH
  if(vdj.gene == 'IGH'){
    class.switch.genes <- c('IGHA2','IGHE','IGHG4','IGHG2','IGHA1',
                            'IGHG1','IGHG3','IGHM')
    output.cs.segments <- output.segments %>%
      dplyr::filter(segment %in% class.switch.genes) %>%
      dplyr::mutate(segment.fraction = !!fraction_col - dplyr::lag(!!fraction_col, default = 0)) %>%
      dplyr::mutate(segment.fraction = ifelse(segment.fraction < 0, 0, segment.fraction)) %>%
      dplyr::mutate(segment = gsub(paste0('_',focal_segment_split[2]),'',segment)) %>%
      dplyr::select(sample, segment, segment.fraction)

    output.all.segments <- rbind(output.cs.segments,output.all.segments)

  }

  # Add in Shannon divergence scores for TRAV and TRAJ segments
  total.fraction <- output.summary %>% dplyr::select(!!fraction_col) %>% `[[`(1)

  if(total.fraction > 0){
    shannon.trav <- output.v.segments %>%
      dplyr::mutate(p = segment.fraction/total.fraction) %>%
      dplyr::filter(p > 0) %>%
      dplyr::mutate(shannon.value = -(p*log(p))) %>%
      dplyr::group_by(sample) %>%
      dplyr::summarise(shannon.trav = sum(shannon.value)) %>%
      dplyr::select(shannon.trav) %>% `[[`(1)

    shannon.traj <- output.j.segments %>%
      dplyr::mutate(p = segment.fraction/total.fraction) %>%
      dplyr::filter(p > 0) %>%
      dplyr::mutate(shannon.value = -(p*log(p))) %>%
      dplyr::group_by(sample) %>%
      dplyr::summarise(shannon.trav = sum(shannon.value)) %>%
      dplyr::select(shannon.trav) %>% `[[`(1)
    if(length(shannon.traj) == 0){
      shannon.traj <- NA}
    if(length(shannon.trav) == 0){
      shannon.trav <- NA}
  }else{
    shannon.trav <- NA
    shannon.traj <- NA
  }


  output.summary[[paste0('shannon.',vseg, '.div')]] <- shannon.trav
  output.summary[[paste0('shannon.',jseg, '.div')]] <- shannon.traj

  # Add in class switch metrics for IGH
  if(vdj.gene == 'IGH'){
    cs.fraction <- sum(output.cs.segments$segment.fraction)
    igA.frac <- output.cs.segments %>%
      dplyr::filter(segment %in% c('IGHA2','IGHA1')) %>%
      dplyr::select(segment.fraction) %>% `[[`(1) %>% sum()
    igG.frac <- output.cs.segments %>%
      dplyr::filter(segment %in% c('IGHG1','IGHG2','IGHG3','IGHG4')) %>%
      dplyr::select(segment.fraction) %>% `[[`(1) %>% sum()
    igE.frac<- output.cs.segments %>%
      dplyr::filter(segment %in% c('IGHE')) %>%
      dplyr::select(segment.fraction) %>% `[[`(1) %>% sum()
    igMD.frac <- output.summary$IGH.bcell.fraction - (igA.frac + igG.frac + igE.frac)
    output.summary[['class.switch.frac']] <- cs.fraction
    output.summary[['igA.frac']] <- igA.frac
    output.summary[['igG.frac']] <- igG.frac
    output.summary[['igE.frac']] <- igE.frac
    output.summary[['igMD.frac']] <- igMD.frac

  }

  return(list(output.summary, output.all.segments, segmodel.out))
}
