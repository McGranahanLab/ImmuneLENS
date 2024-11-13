#' Function to run ImmuneLENS for WGS for individual V and J segments
#'
#' @param test.logR data frame containing logR values to be fitted to constrained model
#' @param vdj.gene V(D)J gene to use
#' @param sample_name name of the sample
#' @param  hg19_or_38 version of the genome
#' @param GC.correct whether GC.correct is being used
#' @param exons Location of exons
#' @param exonList GC content of exons
#' @param gene.fasta FASTA file for VDJ gene e.g TCRA
#' @param gene.fasta.start offset for FASTA file, e.g 21999999 for TCRA in hg19
#' @param sliding number of bp for gc windows (default = 1000)
#' @param IGH.gc.constraint TRUE/FALSE, use additional constraints for GC correction smooth.gc and smooth.gc2 in the IGH solution (default TRUE to prevent over-fitting)
#' @param IGH.gc.constraint.value Value to constrain smooth.gc and smooth.gc2 by (default 0.01) 
#' @param TCRB.gc.constraint TRUE/FALSE, use additional constraints for GC correction smooth.gc and smooth.gc2 in the TCRB solution (default TRUE to prevent over-fitting)
#' @param allGC TRUE/FALSE, when TRUE make model where everything is due to GC
#' @return data frame of TCRA T cell fractions for VDJ segments
#' @name getVDJfraction_ImmuneLENS
#'
getVDJfraction_ImmuneLENS <- function(test.logR, vdj.gene, sample_name,
                                        hg19_or_38 = 'hg38', GC.correct = TRUE,
                                        exons = NULL,
                                        exonList = NULL,
                                        gene.fasta = NULL,
                                        gene.fasta.start = NULL,
                                        sliding = 1000,
                                        IGH.gc.constraint = TRUE,
                                        IGH.gc.constraint.value = 0.01,
                                        TCRB.gc.constraint = TRUE,
                                        allGC = FALSE){

  # Solve binding issues
  end <- start <- segName <- NULL
  Ratio <- exon.gc2 <- smooth.gc2 <- comparison <- NULL
  segName2 <- constraint <- logR <- segment <- VDJ.tcell.fraction <- NULL
  X2 <- NULL
  pos <- GC <- exon.gc <- smooth.gc <- NULL
  
  # Set some values:
  absval_val = 1e-5
  
  exclude.segs <- c('TRDD1','TRDD2','TRDD3','TRDC','TRAV11','TRAV8_4','TRAV8_5',
                    'TRAV15','TRAV31','TRAV7','TRAV9_1','TRAV18','TRAV8_7',
                    'TRAV28','TRAV32','TRAV33','TRAV37','TRDV2','TRDJ1','TRDJ2',
                    'TRDJ3','TRDJ4','TRAJ55','TRAJ51','TRAJ2','TRAJ1','TRAJ60',
                    'TRAJ59','TRAJ24','TRAJ14','TRBV7_1','TRBV6_2','TRBV8_1',
                    'TRBV5_2','TRBV8_2','TRBV12_3','TRBVA','TRBV26','TRBVB',
                    'TRBV1','TRBV5_3','TRB10_1','TRBV7_4','TRBV7_5','TRBV6_7',
                    'TRBV7_7','TRBV6_8','TRBV5_7','TRBV17','TRBV22_1','TRBC2',
                    'TRBV30','TRBD1','TRBC1','TRGC2','TRGC1',
                    'IGHGP','IGHEP1','IGHD',
                    'IGHD1_26', 'IGHD6_25', 'IGHD5_24','IGHD4_23',
                    'IGHD3_22','IGHD2_21', 'IGHD1_20', 'IGHD6_19',
                    'IGHD5_18', 'IGHD4_17','IGHD3_16', 'IGHD2_15', 'IGHD1_14',
                    'IGHD6_13', 'IGHD5_12','IGHD4_11', 'IGHD3_10', 'IGHD3_9',
                    'IGHD2_8', 'IGHD1_7','IGHD6_6' , 'IGHD5_5', 'IGHD4_4',
                    'IGHD3_3', 'IGHD2_2', 'IGHD1_1','IGHD7_27',
                    # IGHV psedogenes:
                    'IGHVII_1_1', 'IGHVIII_2_1', 'IGHVIII_5_1', 'IGHVIII_5_2', 
                    'IGHVIII_11_1', 'IGHVIII_13_1', 'IGHVII_15_1', 'IGHVIII_16_1',
                    'IGHVII_22_1', 'IGHVIII_22_2', 'IGHVIII_25_1', 'IGHVIII_26_1',
                    'IGHVII_26_2', 'IGHVII_28_1', 'IGHVII_30_1', 'IGHVII_30_21',
                    'IGHVII_33_1', 'IGHVIII_38_1', 'IGHVII_40_1', 'IGHVII_43_1',
                    'IGHVIII_44', 'IGHVIV_44_1', 'IGHVII_44_2', 'IGHVII_46_1',
                    'IGHVIII_47_1', 'IGHVII_49_1', 'IGHVII_51_2', 'IGHVII_53_1',
                    'IGHVII_60_1', 'IGHVII_62_1', 'IGHVII_65_1', 'IGHVII_67_1', 
                    'IGHVIII_67_2', 'IGHVIII_67_3', 'IGHVIII_67_4', 'IGHVII_74_1',
                    'IGHVIII_76_1', 'IGHVII_78_1', 'IGHVIII_82', 'IGHVII_20_1',
                    'IGHVII_31_1', 'IGHVIII_51_1',
                    'IGLC7','IGLL5',
                    'IGLC1','IGLC2','IGLC3','IGLC4','IGLC5','IGLC6','IGKC')
  
  class.switch.genes <- c('IGHA2','IGHE','IGHG4','IGHG2','IGHA1',
                          'IGHG1','IGHG3','IGHM')
  
  norm_col_names <- character()
  
  if(GC.correct){
    # Calculate the exon GC values (note this could be saved internally + loaded to save)
    TCRA.gc.df <- slidingwindowplot_alt(sliding, gene.fasta[[1]])
    TCRA.gc.df$pos <- TCRA.gc.df$loc + gene.fasta.start
    gam.model <- mgcv::gam(GC~s(pos, bs = 'cs'), data = TCRA.gc.df)

    # This is the function to get the smoothed GC content at any position
    get_gc_prediction <- function(x){
      mgcv::predict.gam(gam.model,newdata = data.frame(pos = x))}

    test.logR <- test.logR %>%
      dplyr::mutate(exon = exonPosFun_v(pos, exons)) %>%
      dplyr::left_join(exonList, 'exon') %>%
      dplyr::rename(exon.gc = GC) %>%
      dplyr::mutate(exon.gc2 = exon.gc^2) %>%
      dplyr::mutate(smooth.gc = get_gc_prediction(pos))  %>%
      dplyr::mutate(smooth.gc2 = smooth.gc^2)

    norm_col_names <- c('exon.gc','exon.gc2','smooth.gc','smooth.gc2')

  }

  segment.ranges.list <- calculateSegmentRanges(vdj.gene = vdj.gene,
                                           hg19_or_38 = hg19_or_38, 
                                           exclude.segs = exclude.segs,
                                           class.switch.genes = class.switch.genes)
  segment.ranges <- segment.ranges.list[[1]]
  new_seg_name <-   segment.ranges.list[[2]]
  
  if(GC.correct){
      test.logR2 <-   test.logR %>%
        dplyr::select(pos, Ratio, exon.gc, exon.gc2,smooth.gc, smooth.gc2)
  }else{
    test.logR2 <-   test.logR %>%
      dplyr::select(pos, Ratio)
  }

  # Make sure all segments have values
  segment_ranges_dt <- data.table::as.data.table(segment.ranges)
  test_logR2_dt <- data.table::as.data.table(test.logR2)
  
  # Rename columns in test_logR2_dt to be compatible with foverlaps
  data.table::setnames(test_logR2_dt, old = "pos", new = "start")
  test_logR2_dt[, end := start]  # Create an end column (same as start, as pos is a single point)
  
  # Set the key for both data.tables
  data.table::setkey(segment_ranges_dt, start, end)
  data.table::setkey(test_logR2_dt, start, end)
  
  # Perform the overlap join
  overlaps <- data.table::foverlaps(test_logR2_dt, segment_ranges_dt, type = "within", nomatch = 0L)
  # segment.positions.lengths <- overlaps[, .N, by = segName]$N
  
  segment.positions.counts <- overlaps[, .N, by = segName]
  segment.ranges.tmp <- segment.ranges %>% 
    dplyr::left_join(segment.positions.counts, 'segName') %>%
    dplyr::mutate(N = ifelse(is.na(N),0,N))
  
  segment.positions.lengths <- segment.ranges.tmp$N
  
  # Update the segment.ranges to remove the ones with 0 positions
  segment.ranges.upd <- segment.ranges
  if(length(which(segment.positions.lengths == 0)) > 0){
    for(i in rev(which(segment.positions.lengths == 0))){
      upd.end <- segment.ranges.upd$end[i]
      segment.ranges.upd <- segment.ranges.upd[-i, ]
      segment.ranges.upd$end[i - 1] <- upd.end
    }
    segment.positions.lengths2 <- segment.positions.lengths[-which(segment.positions.lengths == 0)]
  }else{
    segment.positions.lengths2 <- segment.positions.lengths
  }

  # For IGH check if the VJ segment has been removed and re-create:
  if(vdj.gene == 'IGH'){
    if(segment.positions.lengths[which(segment.ranges$segName == new_seg_name)] == 0){
      remaining_segs <- segment.ranges.upd$segName
      segment.ranges.upd <- segment.ranges.upd %>%
        dplyr::mutate(segType = ifelse(grepl('V',segName),'V','J')) %>%
        dplyr::mutate(next_segType = ifelse(grepl('V',dplyr::lead(segName)),'V','J')) %>%
        dplyr::mutate(prev_segType = ifelse(grepl('V',dplyr::lag(segName)),'V','J')) 
      
      v_seg_loc <- which(segment.ranges.upd$next_segType == 'V')
      new_seg_name <- paste0(segment.ranges.upd$segName[v_seg_loc[1]],'_',
                             segment.ranges.upd$segName[v_seg_loc[1] + 1])
      
      segment.ranges.upd$segName[v_seg_loc[1]] <- new_seg_name
      segment.ranges.upd  <- segment.ranges.upd [-c(v_seg_loc[1] + 1), ]
      segment.ranges.upd  <- segment.ranges.upd  %>%
        dplyr::select(segName, start, end) %>%
        dplyr::mutate(segName = as.character(segName))
    }
  } 
  
  # Add columns to the matrix
  test_logR2_dt <- data.table::as.data.table(test.logR2)
  segment_ranges_dt <- data.table::as.data.table(segment.ranges.upd)
  norm_cols <- setdiff(colnames(test.logR2),c('pos','Ratio'))
  
  if(length(norm_cols) == 0){
    merged_data <- test_logR2_dt[, .(pos,Ratio)]
  }else{
    merged_data <- test_logR2_dt[, c('pos','Ratio', ..norm_cols)]
  }
  
  # Initialize segment columns with 0
  for (seg in unique(segment_ranges_dt$segName)) {
    merged_data[, (seg) := 0L]
  }
  
  # For each segment, set the indicator to 1 for positions within the segment
  data.table::setnames(test_logR2_dt, old = "pos", new = "start")
  test_logR2_dt[, end := start]  # Create an end column (same as start, as pos is a single point)
  
  # Set the key for both data.tables
  
  # Prepare for foverlaps
  data.table::setkey(test_logR2_dt, start, end)
  data.table::setkey(segment_ranges_dt, start, end)
  
  # Perform overlap join
  overlaps <- data.table::foverlaps(test_logR2_dt, segment_ranges_dt, type = "within", nomatch = NULL)
  data.table::setnames(overlaps, old = "i.start", new = "pos")
  
  # Keep relevant columns
  overlaps <- overlaps[, .(pos, segName)]
  
  positions_list <- split(overlaps$pos, overlaps$segName)
  for (seg in unique(segment_ranges_dt$segName)) {
    # Set the indicator to 1 for these positions
    merged_data[pos %in% positions_list[[seg]], (seg) := 1L]
  }
  
  # Convert back to data.frame if needed
  test.logR2 <- as.data.frame(merged_data)
  
  colnames(test.logR2) <- gsub('-','_',colnames(test.logR2))
  segment.ranges.upd$segName <- gsub('-','_',segment.ranges.upd$segName)

  # Remove IGHM column if there
  segment.ranges.names <- segment.ranges.upd$segName

  if(vdj.gene == 'IGH'){
    ighm.col <- which(colnames(test.logR2) == 'IGHM')
    if(length(ighm.col) > 0){
      test.logR2 <- test.logR2[,-ighm.col]
    }
    segment.ranges.names <- setdiff(segment.ranges.names, 'IGHM')
  }

  rm(list = c('merged_data', 'overlaps', 'positions_list',
              'segment_ranges_dt', 'test_logR2_dt'))
  gc()
  
  if(GC.correct){
      model.lm <- paste0('Ratio ~ 0 + exon.gc + exon.gc2 + smooth.gc + smooth.gc2 + ',
                         (paste0(segment.ranges.names, collapse = ' + ')))
  }else{
    model.lm <- paste0('Ratio ~ 0 + ',
                       (paste0(segment.ranges.names, collapse = ' + ')))
  }

  test.lm.unconstrained <- lm(model.lm, test.logR2[,-1])
  
  myConstraints <- generateConstraints(segment.ranges.upd = segment.ranges.upd,
                                       vdj.gene = vdj.gene,
                                       GC.correct = GC.correct,
                                       class.switch.genes = class.switch.genes,
                                       IGH.gc.constraint = IGH.gc.constraint,
                                       IGH.gc.constraint.value = IGH.gc.constraint.value ,
                                       TCRB.gc.constraint = TCRB.gc.constraint,
                                       allGC = allGC)

  # Prepare for custom fitting functions
  y <- as.matrix(test.lm.unconstrained$model[, attr(test.lm.unconstrained$terms, "response")])
  so <- summary(test.lm.unconstrained)
  weights <- weights(test.lm.unconstrained)
  p <- length(stats::coef(test.lm.unconstrained))
  b.unrestr <- stats::coef(test.lm.unconstrained)
  Sigma <- stats::vcov(test.lm.unconstrained)
  residuals <- test.lm.unconstrained$residuals
  
  constraint.OUT <- get_constraint_mat(test.lm.unconstrained, constraints =  myConstraints)
  
  terms <- test.lm.unconstrained$terms
  fitted <- test.lm.unconstrained$fitted
  df.residual <- test.lm.unconstrained$df.residual
  rm(test.lm.unconstrained)
  gc()
  
  X <- as.matrix(test.logR2[,-c(1,2)])
  
  restr.lm <- conLM.lm_custom(y = y, X = X, so = so, weights = weights,p = p,
                          segment.positions.lengths = segment.positions.lengths2,
                          b.unrestr = b.unrestr, Sigma = Sigma,
                          residuals = residuals, constraint.OUT = constraint.OUT,
                          terms = terms, fitted = fitted, df.residual = df.residual,
                          control = list(absval = absval_val),
                          norm_cols = norm_col_names)

  num.par <- length(restr.lm$parTable$est)
  num.obs <- dim(test.logR2)[1]

  vdj_segment_output <- data.frame(logR = restr.lm$b.restr) %>%
    tibble::rownames_to_column('segment') %>%
    dplyr::mutate(VDJ.tcell.fraction = 1- 2^logR) %>%
    dplyr::mutate(sample = sample_name) %>%
    dplyr::select(sample, segment, logR, VDJ.tcell.fraction) %>%
    dplyr::mutate(modelR2 =  restr.lm$R2.reduced) %>%
    dplyr::mutate(model.loglik =  restr.lm$loglik) %>%
    dplyr::mutate(model.s2=  restr.lm$s2) %>%
    dplyr::mutate(model.BIC = -2*restr.lm$loglik + log(num.obs)*num.par)

  colnames(vdj_segment_output)[4] <- ifelse(vdj.gene %in% c('TCRA','TCRB','TCRG'),
                                            paste0(vdj.gene,'.tcell.fraction'),
                                            paste0(vdj.gene,'.bcell.fraction'))

  # Remove some of the large objects
  rm(restr.lm)
  gc()

  return(vdj_segment_output)
}


generateConstraints <- function(segment.ranges.upd, vdj.gene, GC.correct = FALSE,
                                allGC = FALSE, IGH.gc.constraint = FALSE, 
                                IGH.gc.constraint.value = 0.00001, TCRB.gc.constraint = FALSE,
                                norm_col_names = NULL, sample_name = NULL,
                                class.switch.genes = NULL){
  # Solve binding issues
  segName <- comparison <- segName2 <- constraint <- NULL
  
   if(vdj.gene == 'IGH'){
    start_segment <- setdiff(segment.ranges.upd$segName, class.switch.genes)[1]
  }else{
    start_segment <- segment.ranges.upd$segName[1]
  }
  
  if(grepl('J',start_segment)){
    myConstraints1 <- segment.ranges.upd %>%
      dplyr::mutate(segName2 = paste0(dplyr::lag(segName,default = '0'),'\n')) %>%
      dplyr::mutate(comparison = ifelse(grepl('J',segName), ' < ',' > ')) %>%
      dplyr::mutate(constraint = paste0(segName, comparison, segName2))    %>%
      dplyr::select(constraint) %>% `[[`(1)
  }else{
    myConstraints1 <- segment.ranges.upd %>%
      dplyr::mutate(segName2 = paste0(dplyr::lag(segName,default = '0'),'\n')) %>%
      dplyr::mutate(comparison = ifelse(grepl('J',segName),
                                        ifelse(grepl('V',segName),' < ',' > '),' < ')) %>%
      dplyr::mutate(constraint = paste0(segName, comparison, segName2))    %>%
      dplyr::select(constraint) %>% `[[`(1)
    
  }
  
  if(vdj.gene == 'IGH'){
    # Fix Class switching comparisons (remove IGHM):
    # Class switch constraint loc:
    class.switch.constraint.loc <- lapply(setdiff(class.switch.genes,'IGHM'),
                                          function(x) grep(x,myConstraints1)) %>% unlist() %>% unique()
    IGHM.constraint.loc <- grep('IGHM',myConstraints1)
    class.switch.constraint.loc <- setdiff(class.switch.constraint.loc, IGHM.constraint.loc)
    if(length(class.switch.constraint.loc) == 0 | length(IGHM.constraint.loc) != 2){
      warning('Not enough bases with coverage in IGH class switch region')
      return(NULL)
    }
    
    
    myConstraints1[class.switch.constraint.loc] <- gsub('>','<',myConstraints1[class.switch.constraint.loc])
    
    final_cs_seg <- segment.ranges.upd$segName[max(which(segment.ranges.upd$segName %in% setdiff(class.switch.genes, 'IGHM')))]
    first_v_seg <- segment.ranges.upd$segName[min(grep('J',segment.ranges.upd$segName))]
    JV_seg <- segment.ranges.upd$segName[min(grep('V',segment.ranges.upd$segName))]
    next_seg_name  <- JV_seg
    myConstraints1[IGHM.constraint.loc[1]] <- paste0(final_cs_seg, ' > ', JV_seg, '\n')
    myConstraints1[IGHM.constraint.loc[2]] <-  paste0(first_v_seg, ' < 0\n')
  }
  
  
  
  myConstraints1 <- myConstraints1 %>% paste0(collapse = ' ')
  
  myConstraints <- paste0('\n ',myConstraints1, ' ',
                          segment.ranges.upd$segName[length(segment.ranges.upd$segName)],
                          ' < 0')
  
  if(vdj.gene == 'IGH' & GC.correct & IGH.gc.constraint){
    # Add in constraints for smooth.gc and smooth.gc2 to prevent over-fitting of IGH solution
    # Try preventing all the fraction from occuring at the very end:
    # IGHV3_73 + IGHV3_74 + IGHV3_75 + IGHV3_76 + IGHV5_78 + IGHV3_79 + IGHV4_80 + IGHV7_81
    if('IGHV3_73' %in% segment.ranges.upd$segName){
      gc.constraints <- paste0('\n IGHV3_73  > ',IGH.gc.constraint.value,'*',JV_seg)
    }else{
      end_v_seg <- utils::tail(segment.ranges.upd$segName,n = 7)[1]
      gc.constraints <- paste0('\n ',end_v_seg,'  > ',IGH.gc.constraint.value,'*',JV_seg)
    }
    
    
    myConstraints <- paste0(myConstraints, gc.constraints)
  }
  
  if(vdj.gene == 'TCRB' & GC.correct & TCRB.gc.constraint){
    # 1) smooth.gc + smooth.gc2 > 0, 2) exon.gc2 + smooth.gc2 > 2
    # 3) abs(ALL) > 3 and < 10
    gc.constraints <- c('\n smooth.gc + smooth.gc2 > 0 \n exon.gc2 + smooth.gc2 > 2 \n exon.gc2 + smooth.gc - exon.gc - smooth.gc2 < 10')
    myConstraints <- paste0(myConstraints, gc.constraints)
  }
  
  
  if((allGC)){
    next_seg_name <- gsub('-','_',next_seg_name)
    new.constraint <- paste0('\n ',next_seg_name, ' > -0.00001')
    myConstraints <- paste0(myConstraints, new.constraint)
  }
  
  return(myConstraints)
  
}


calculateSegmentRanges <- function(vdj.gene, hg19_or_38, 
                                  exclude.segs,
                                  class.switch.genes){
  # Solve binding issues
  segName <- segType <- end <- prev_segType <- NULL
  start <- next_segType <- start2 <- end2 <-NULL
  hgnc_symbol <- start_position_hg38 <- end_position_hg38 <- NULL
  start_position_hg19 <- end_position_hg19 <- NULL
 
  new_seg_name <- NA
  
  if(vdj.gene == 'TCRA'){
    if(hg19_or_38 == 'hg38'){
      segment.ranges <- tcra_seg_hg38_vdj[-c(1:4),]
    }else{
      if(hg19_or_38 == 'hg19'){
        # data("tcra_seg_hg19_vdj_version")
        segment.ranges <- tcra_seg_hg19_vdj[-c(1:4),]
      }else{
        stop('hg19_or_38 must be hg19 or hg38')
      }
    }
    segment.ranges <- segment.ranges %>%
      dplyr::mutate(segName = gsub('-','_', segName)) %>%
      dplyr::filter(!segName %in% exclude.segs) %>%
      dplyr::mutate(segType = ifelse(grepl('V',segName),'V','J')) %>%
      dplyr::mutate(next_segType = ifelse(grepl('V',dplyr::lead(segName)),'V','J')) %>%
      dplyr::mutate(prev_segType = ifelse(grepl('V',dplyr::lag(segName)),'V','J')) %>%
      dplyr::mutate(start2 = ifelse(segType == 'V',end + 1,
                                    ifelse(prev_segType == 'V',
                                           dplyr::lag(end) + 1,dplyr::lag(start)))) %>%
      dplyr::mutate(end2 = ifelse(segType == 'V',
                                  ifelse(next_segType == 'J',
                                         dplyr::lead(start) - 1,dplyr::lead(end)),
                                  start -1))
    
    j_seg_loc <- which(segment.ranges$next_segType == 'J')
    new_seg_name <- paste0(segment.ranges$segName[j_seg_loc[1]],
                           '_', segment.ranges$segName[j_seg_loc[1] + 1])
    segment.ranges$segName[j_seg_loc[1]] <- new_seg_name
    segment.ranges <- segment.ranges[-c(j_seg_loc[1] + 1), ]
    segment.ranges <- segment.ranges %>%
      dplyr::select(segName, start = start2, end = end2) %>%
      dplyr::mutate(segName = as.character(segName))
    
  }else{
    if(hg19_or_38 == 'hg38'){
      if(vdj.gene %in% c('TCRG','IGH')){
        segment.ranges <- vdj.segments.list[[vdj.gene]] %>%
          dplyr::select(segName = hgnc_symbol, start = start_position_hg38,
                        end = end_position_hg38) %>%
          dplyr::mutate(segName = gsub('-','_', segName)) %>%
          dplyr::filter(!is.na(start)) %>%
          dplyr::filter(!segName %in% exclude.segs) %>%
          dplyr::mutate(segType = ifelse(grepl('V',segName),'V','J')) %>%
          dplyr::mutate(next_segType = ifelse(grepl('V',dplyr::lead(segName)),'V','J')) %>%
          dplyr::mutate(prev_segType = ifelse(grepl('V',dplyr::lag(segName)),'V','J')) %>%
          dplyr::mutate(start2 = ifelse(segType == 'J',end + 1,
                                        ifelse(prev_segType == 'J',
                                               dplyr::lag(end) + 1,dplyr::lag(start)))) %>%
          dplyr::mutate(end2 = ifelse(segType == 'J',
                                      ifelse(next_segType == 'V',
                                             dplyr::lead(start) - 1,dplyr::lead(end)),
                                      start -1))
        
        if(vdj.gene == 'IGH'){
          # Define class switching gene segments and update segment ranges
          
          segment.ranges <- segment.ranges %>%
            dplyr::mutate(segType = ifelse(segName %in% class.switch.genes,
                                           'CS',segType)) %>%
            dplyr::mutate(next_segType= ifelse(segName %in% class.switch.genes[-8],
                                               'CS',next_segType)) %>%
            dplyr::mutate(prev_segType = ifelse(segName %in% c(class.switch.genes,'IGHJ6'),
                                                'CS',next_segType))
          
          # Class switch segment guide:
          # Segments IGHA2 to IGHG3 monotonically decreasing
          # IGHG3 segment can not be lower than minimal fraction
          # (e.g. class switched B cells can not be greater than all B cells)
          # IGHM segment back to 0
          # VDJ segments then same as before
          
        }
        
        
        v_seg_loc <- which(segment.ranges$next_segType == 'V')
        new_seg_name <- paste0(segment.ranges$segName[v_seg_loc[1]],'_',
                               segment.ranges$segName[v_seg_loc[1] + 1])
        segment.ranges$segName[v_seg_loc[1]] <- new_seg_name
        segment.ranges <- segment.ranges[-c(v_seg_loc[1] + 1), ]
        segment.ranges <- segment.ranges %>%
          dplyr::select(segName, start = start2, end = end2) %>%
          dplyr::mutate(segName = as.character(segName))
      }else{
        segment.ranges <- vdj.segments.list[[vdj.gene]] %>%
          dplyr::select(segName = hgnc_symbol, start = start_position_hg38,
                        end = end_position_hg38) %>%
          dplyr::mutate(segName = gsub('-','_', segName)) %>%
          dplyr::filter(!is.na(start)) %>%
          dplyr::filter(!segName %in% exclude.segs) %>%
          dplyr::mutate(segType = ifelse(grepl('V',segName),'V','J')) %>%
          dplyr::mutate(next_segType = ifelse(grepl('V',dplyr::lead(segName)),'V','J')) %>%
          dplyr::mutate(prev_segType = ifelse(grepl('V',dplyr::lag(segName)),'V','J')) %>%
          dplyr::mutate(start2 = ifelse(segType == 'V',end + 1,
                                        ifelse(prev_segType == 'V',
                                               dplyr::lag(end) + 1,dplyr::lag(start)))) %>%
          dplyr::mutate(end2 = ifelse(segType == 'V',
                                      ifelse(next_segType == 'J',
                                             dplyr::lead(start) - 1,dplyr::lead(end)),
                                      start -1))
        
        J_seg_loc <- which(segment.ranges$next_segType == 'J')
        new_seg_name <- paste0(segment.ranges$segName[J_seg_loc[1]],'_',
                                segment.ranges$segName[J_seg_loc[1] + 1])
        segment.ranges$segName[J_seg_loc[1]] <- new_seg_name
        segment.ranges <- segment.ranges[-c(J_seg_loc[1] + 1), ]
        segment.ranges <- segment.ranges %>%
          dplyr::select(segName, start = start2, end = end2) %>%
          dplyr::mutate(segName = as.character(segName))
      }
    }else{
      if(hg19_or_38 == 'hg19'){
        if(vdj.gene %in% c('TCRG','IGH')){
          segment.ranges <- vdj.segments.list[[vdj.gene]] %>%
            dplyr::select(segName = hgnc_symbol, start = start_position_hg19,
                          end = end_position_hg19) %>%
            dplyr::mutate(segName = gsub('-','_', segName)) %>%
            dplyr::filter(!is.na(start)) %>%
            dplyr::filter(!segName %in% exclude.segs) %>%
            dplyr::mutate(segType = ifelse(grepl('V',segName),'V','J')) %>%
            dplyr::arrange(start) %>%
            dplyr::mutate(next_segType = ifelse(grepl('V',dplyr::lead(segName)),'V','J')) %>%
            dplyr::mutate(prev_segType = ifelse(grepl('V',dplyr::lag(segName)),'V','J')) %>%
            dplyr::mutate(start2 = ifelse(segType == 'J',end + 1,
                                          ifelse(prev_segType == 'J',
                                                 dplyr::lag(end) + 1,dplyr::lag(start)))) %>%
            dplyr::mutate(end2 = ifelse(segType == 'J',
                                        ifelse(next_segType == 'V',
                                               dplyr::lead(start) - 1,dplyr::lead(end)),
                                        start -1))
          
          V_seg_loc <- which(segment.ranges$next_segType == 'V')
          new_seg_name <- paste0(segment.ranges$segName[V_seg_loc[1]],'_',
                                  segment.ranges$segName[V_seg_loc[1] + 1])
          
          segment.ranges$segName[V_seg_loc[1]] <- new_seg_name
          segment.ranges <- segment.ranges[-c(V_seg_loc[1] + 1), ]
          segment.ranges <- segment.ranges %>%
            dplyr::select(segName, start = start2, end = end2) %>%
            dplyr::mutate(segName = as.character(segName))
        }else{
          segment.ranges <- vdj.segments.list[[vdj.gene]] %>%
            dplyr::select(segName = hgnc_symbol, start = start_position_hg19,
                          end = end_position_hg19) %>%
            dplyr::mutate(segName = gsub('-','_', segName)) %>%
            dplyr::filter(!is.na(start)) %>%
            dplyr::filter(!segName %in% exclude.segs) %>%
            dplyr::mutate(segType = ifelse(grepl('V',segName),'V','J')) %>%
            dplyr::arrange(start) %>%
            dplyr::mutate(next_segType = ifelse(grepl('V',dplyr::lead(segName)),'V','J')) %>%
            dplyr::mutate(prev_segType = ifelse(grepl('V',dplyr::lag(segName)),'V','J')) %>%
            dplyr::mutate(start2 = ifelse(segType == 'V',end + 1,
                                          ifelse(prev_segType == 'V',
                                                 dplyr::lag(end) + 1,dplyr::lag(start)))) %>%
            dplyr::mutate(end2 = ifelse(segType == 'V',
                                        ifelse(next_segType == 'J',
                                               dplyr::lead(start) - 1,dplyr::lead(end)),
                                        start -1))
          
          J_seg_loc <- which(segment.ranges$next_segType == 'J')
          new_seg_name <- paste0(segment.ranges$segName[J_seg_loc[1]],'_',
                                  segment.ranges$segName[J_seg_loc[1] + 1])
          segment.ranges$segName[J_seg_loc[1]] <- new_seg_name
          segment.ranges <- segment.ranges[-c(J_seg_loc[1] + 1), ]
          segment.ranges <- segment.ranges %>%
            dplyr::select(segName, start = start2, end = end2) %>%
            dplyr::mutate(segName = as.character(segName))
        }
        
      }
    }
  }
  return(list(segment.ranges, new_seg_name))
  
}

conLM.lm_custom <- function (y, X, so, weights,p,b.unrestr,Sigma,residuals,
                             constraint.OUT,terms,
                             fitted, df.residual,
                             norm_cols,  control = list()) 
{
  Amat <- NULL
  bvec <- NULL
  meq <- 0L
  s2 <- so$sigma^2
  b.unrestr[abs(b.unrestr) < ifelse(is.null(control$tol), sqrt(.Machine$double.eps), 
                                    control$tol)] <- 0L
  n <- dim(X)[1]
  object.restr <- list(residuals = residuals, weights = weights)

  parTable <- constraint.OUT$partable 

  Amat <- constraint.OUT$Amat # Different - needs to be calculated

  out.solver <- con_solver_lm_custom(X = X, y = y, w = weights, 
                                   Amat = Amat, meq = meq,
                                   norm_cols = norm_cols,
                                   absval = ifelse(is.null(control$absval),
                                                   sqrt(.Machine$double.eps), control$absval),
                                   maxit = ifelse(is.null(control$maxit), 10000, control$maxit),
                                   b.unrestr = b.unrestr)
    
    b.restr <- out.solver$qp$solution
    names(b.restr) <- names(b.unrestr)
    b.restr[abs(b.restr) < ifelse(is.null(control$tol), sqrt(.Machine$double.eps), 
                                  control$tol)] <- 0L
    
      fitted <- X %*% b.restr
      residuals <- y - fitted
      object.restr <- list(residuals = residuals, weights = weights)
      # Rewrite:
      ll.restr <- get_lm_loglik(object.restr)
      
      if (is.null(weights)) {
        mss <- if (attr(terms, "intercept")) {
          sum((fitted - mean(fitted))^2)
        }else {
          sum(fitted^2)
        }
        rss <- sum(residuals^2)
      }else {
        mss <- if (attr(terms, "intercept")) {
          m <- sum(weights * fitted/sum(weights))
          sum(weights * (fitted - m)^2)
        }
        else {
          sum(weights * fitted^2)
        }
        rss <- sum(weights * residuals^2)
      }
      R2.reduced <- mss/(mss + rss)
      if (is.null(weights)) {
        s2 <- sum(residuals^2)/df.residual
      } else {
        s2 <- sum(weights * residuals^2)/df.residual
      }

    OUT <- list( parTable = parTable, b.restr = b.restr, 
                 R2.reduced = so$r.squared, s2 = s2, loglik = ll.restr)
  return(OUT)
}

# Custom solver for linear regression with constraints (fast):
con_solver_lm_custom <- function (X, y, w = NULL, Amat, meq, maxit = 10000,
                              absval = sqrt(.Machine$double.eps),
                              b.unrestr,
                              norm_cols) 
{
  bvec <- rep(0, nrow(Amat))
  val <- 0
  X <- as.matrix(X)
  y <- as.matrix(y)
  # b.restr <- c(coef(lm.fit(x = X, y)))
  b.restr <- b.unrestr
  norm_loc <- which(colnames(X) %in% norm_cols)
  n <- nrow(X[,-norm_loc])
  p <- ncol(X[,-norm_loc])
  tmp_X <- X[,-norm_loc]
  first_entry <- sapply(seq(p), function(x){which(tmp_X[,x] == 1)[1]})
  last_entry <- utils::tail(which(tmp_X[,p] == 1))[1]
  rm(tmp_X)
  gc()
  rep_numbers <- rep(0,p+2)
  for(i in seq(p)){
    if(i == 1) {
      rep_numbers[i] <- first_entry[i] - 1
    }else{
      rep_numbers[i] <- first_entry[i] - first_entry[i-1] 
    }
  }
  rep_numbers[p + 1] <- last_entry - first_entry[p] + 1
  rep_numbers[p + 2] <- n - last_entry
  
  t_X_y <- t(X) %*% y
  t_X_X <- (t(X) %*% X) 
  
  for (i in 1:maxit) {
    y_est <- lapply(seq(p+2), function(x) rep(c(0, b.restr, 0)[x], rep_numbers[x])) %>% unlist()
    resid <- y - y_est
    if (!is.null(w)) {
      s2 <- sum(w * resid^2) / (n - p)
      yVx <- (1 / s2) * t(X) %*% (w * y)
      Dmat <- 2 * (1 / s2) * (t(X) %*% (w * X))
    } else {
      s2 <- sum(resid^2) / (n - p)
      # Can these calculations be made more efficient considering the structure of the matrix?
      yVx <- (1 / s2) * t_X_y 
      Dmat <- 2 * (1 / s2) * t_X_X
    }
    dvec <- 2 * yVx
    out.qp <- quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = t(Amat), 
                                  bvec = bvec, meq = meq)
    b.restr <- out.qp$solution
    if (abs(out.qp$value - val) <= absval) {
      break
    } else {
      val <- out.qp$value
    }
    if (i == maxit & abs(out.qp$value - val) > absval) {
      warning(gettextf("'quadprog' failed to converge in %d steps", 
                       maxit), domain = NA)
    }
  }
  list(qp = out.qp, s2 = s2, Niter = i)
}


get_lm_loglik <- function(lm.out){
  residuals <- as.matrix(lm.out$residuals)
  s2 <- crossprod(residuals)/nrow(residuals)
  n <- length(residuals)
  weights <- rep.int(1, n)

  loglik <- 0.5 * (sum(log(weights)) - n * (log(2 * pi) + 1 - log(n) + 
                                       log(sum(weights * residuals^2))))

  return(loglik)
  
}

get_constraint_mat <- function(model, constraints){
  objectTerms <- stats::terms(model)
  responseIndex <- attr(objectTerms, "response")
  varNames <- as.character(attr(objectTerms, "variables"))[-1]
  responseName <- varNames[responseIndex]
  predCoef <- stats::coef(model)
  predNames <- names(predCoef)
  
  lhs <- rep(responseName, length(predNames))
  op <- rep("~", length(predNames))
  rhs <- predNames
  
  partable <- list(lhs = lhs, op = op, rhs = rhs)
  partable$est <- as.numeric(predCoef)
  partable$label <- predNames
  
  # Check for TCRB/TCRG/IGH constraints
  constraints <- unlist(constraints)
  constraints_usr <- constraints
  constraints  <- gsub(' ','',constraints)
  
  
  CON <- lavaan::lav_constraints_parse(constraints = constraints,
                                       partable = partable, 
                                       theta = partable$est)
  
  CON$constraints <- constraints_usr
  FLAT <- lavaan::lavParseModelString(constraints)
  CON_FLAT <- attr(FLAT, "constraints")
  
  LIST <- list()
  lhs <- unlist(lapply(CON_FLAT, "[[", "lhs"))
  op <- unlist(lapply(CON_FLAT, "[[", "op"))
  rhs <- unlist(lapply(CON_FLAT, "[[", "rhs"))
  LIST$lhs <- lhs
  LIST$op <- op
  LIST$rhs <- c(LIST$rhs, rhs)
  partable$lhs <- c(partable$lhs, LIST$lhs)
  partable$op <- c(partable$op, LIST$op)
  partable$rhs <- c(partable$rhs, LIST$rhs)
  partable$label <- c(partable$label, rep("", length(lhs)))
  
  Amat <- rbind(CON$ceq.JAC, CON$cin.JAC)
   
  OUT <- list(CON = CON, partable = partable, Amat = Amat)
  return(OUT)
}
  

