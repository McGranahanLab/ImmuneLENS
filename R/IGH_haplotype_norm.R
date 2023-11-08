#' Function to normalise IGH coverage based on Haplotype (Experimental!)
#'
#' @param input_cov data frame output from runTcellExtrect_WGS containing segment fractions
#' @param kb_len_threshold If detected germline CNA is longer than this normalise
#' @return List of updated coverage values + regions normalised
#' @importFrom ggplot2 ggplot aes geom_polygon geom_text coord_equal theme_bw theme element_blank
#' @name IGH_haplotype_norm
#' @export
IGH_haplotype_norm <- function(input_cov, kb_len_threshold = 5, ranges_to_normalise = NULL){
 # Get ranges where we see duplications/deletions
  # Issues:
  # 1. Assumes diploid
  # 2. Ignores allelic imbalance
  # 3. Ignores allelic exclusion
  # 4. Assumes B cell content is low (e.g. < 10%) and does not interfere with haplotype calling
  
  # Possible bugs/edge cases to check:
  # empty ranges to normalies
  # High B cell fraction
  if(is.null(ranges_to_normalise)){
    ranges_to_normalise <- calc_haplotype_regions_de_novo_diploid(input_cov, kb_len_threshold)
  }
 
  if(dim(ranges_to_normalise)[1] == 0){
    return(list(input_cov, ranges_to_normalise))
  }
  update_cov <- IGH_haplotype_norm_for_region_diploid(input_cov, ranges_to_normalise)
  
  return(list(update_cov, ranges_to_normalise))
}

IGH_gc_correct <- function(input_cov, hg19_or38 = 'hg38'){
  # Run GC correction on coverage values before running haplotype norm
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
  exons.gc.content <- ImmuneLENS:::exonwindowplot2(IGH.exons.loc, VDJ_fasta[[1]],0)
  
  gc.df <- ImmuneLENS:::slidingwindowplot_alt(1000, VDJ_fasta[[1]])
  gc.df$pos <- gc.df$loc + exon.adjust.loc
  gam.model <- mgcv::gam(GC~s(pos, bs = 'cs'), data = gc.df)
    
  get_gc_prediction <- function(x){
    mgcv::predict.gam(gam.model, newdata = data.frame(pos = x))}
  
  input_cov_gc <- input_cov %>%
    dplyr::mutate(exon = ImmuneLENS:::exonPosFun_v(pos, IGH.exons)) %>%
    dplyr::left_join(exons.gc.content, 'exon') %>%
    dplyr::rename(exon.gc = GC) %>%
    dplyr::mutate(exon.gc2 = exon.gc^2) %>%
    dplyr::mutate(smooth.gc = get_gc_prediction(pos)) %>%
    dplyr::mutate(smooth.gc2 = smooth.gc^2)
  
  # Finally normalise by the GC content
  gc.lm <- lm(reads ~ exon.gc + exon.gc2 + smooth.gc + smooth.gc2,
              y = TRUE, data = input_cov_gc)
  gc.lm.df <- data.frame(pos = input_cov_gc$pos[as.numeric(names(gc.lm$residuals))],
                         reads.gc.correct = gc.lm$residuals)
  if(min(gc.lm.df$reads.gc.correct) < 0){
    gc.lm.df$reads.gc.correct <- gc.lm.df$reads.gc.correct + abs(min(gc.lm.df$reads.gc.correct))
  }

  cov.df.gc.update <- gc.lm.df %>%
    rename(reads = reads.gc.correct)
  # Finally make sure old zeros are still zeros
  zero_pos <- input_cov %>% filter(reads == 0) %>% select(pos) %>% `[[`(1)
  zero_peak <- cov.df.gc.update %>% 
    filter(pos %in% zero_pos) %>% select(reads) %>% `[[`(1) %>% 
    median(na.rm = TRUE)

  cov.df.gc.update$reads <- cov.df.gc.update$reads - zero_peak
  cov.df.gc.update <- cov.df.gc.update %>% 
    mutate(reads = ifelse(reads < 0, 0, reads))
  
  
  return(cov.df.gc.update)
}

IGH_gc_vdj_correct <- function(input_cov, model.fit, hg19_or38 = 'hg38'){
  # Use fit of ImmuneLENS to normalise model
  
  # For input_cov we need 1. LogR values, 2. GC values, 3. segment binary values
  # Get GC:
  # Run GC correction on coverage values before running haplotype norm
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
  exons.gc.content <- ImmuneLENS:::exonwindowplot2(IGH.exons.loc, VDJ_fasta[[1]],0)
  
  gc.df <- ImmuneLENS:::slidingwindowplot_alt(1000, VDJ_fasta[[1]])
  gc.df$pos <- gc.df$loc + exon.adjust.loc
  gam.model <- mgcv::gam(GC~s(pos, bs = 'cs'), data = gc.df)
  
  get_gc_prediction <- function(x){
    mgcv::predict.gam(gam.model, newdata = data.frame(pos = x))}
  
  input_cov_gc <- input_cov %>%
    dplyr::mutate(exon = ImmuneLENS:::exonPosFun_v(pos, IGH.exons)) %>%
    dplyr::left_join(exons.gc.content, 'exon') %>%
    dplyr::rename(exon.gc = GC) %>%
    dplyr::mutate(exon.gc2 = exon.gc^2) %>%
    dplyr::mutate(smooth.gc = get_gc_prediction(pos)) %>%
    dplyr::mutate(smooth.gc2 = smooth.gc^2)
  
  # Get LogR - note do not remove any flagged regions this time
  median.k = 50
  median.thresh = 0
  exons.to.use <- seq_len(dim(IGH.exons)[[1]])
  median.exon.output <- ImmuneLENS:::medianExonCoverage(input_cov, exons.selected,
                                           median.k, median.thresh,
                                           exons.to.use)
  vdj.region.df.filt.exons.median <- median.exon.output[[1]]
  exon.remove <- median.exon.output[[2]]
  
  # Calculate log ratio
  minCov = 0
  col_input <- 'reads'
  col.sym <- rlang::sym(col_input)
  # For random locations with no VDJ effect use beginning and end of TCRA
  tumour.random.covs1 <- vdj.region.df.filt.exons.median  %>%
    dplyr::filter(pos <= vdj.seg[3,]$end) %>%
    dplyr::filter(reads >=  minCov) %>% dplyr::select(!!col.sym) %>% `[[`(1)
  tumour.random.covs2 <- vdj.region.df.filt.exons.median  %>%
    dplyr::filter(pos >= vdj.seg[4,]$start) %>%
    dplyr::filter(reads >=  minCov) %>% dplyr::select(!!col.sym) %>% `[[`(1)
  tumour.random.covs <- c(tumour.random.covs1, tumour.random.covs2)
  # Use median of these values for normalisation to get "logR"
  n1 <- median(tumour.random.covs)
  
  input_cov_gc <-   input_cov_gc %>%
    dplyr::mutate(Ratio = log2(reads/n1))
  
  # Add binary segments from model
  model.fit2 <- model.fit %>%
    dplyr::select(segment, logR) %>%
    dplyr::filter(!segment %in% c('exon.gc','exon.gc2','smooth.gc','smooth.gc2'))
  
  create_col_seg <- function(input_df, start_pos, end_pos,col.name){
    input_df <- input_df %>%
      dplyr::mutate(X1= ifelse(pos >= start_pos & pos <= end_pos, 1, 0))
    colnames(input_df)[which(colnames(input_df) == 'X1')] <- col.name
    return(input_df)
  }
  
  exclude.segs <- c('IGHGP','IGHEP1','IGHD',
                    'IGHD1_26', 'IGHD6_25', 'IGHD5_24','IGHD4_23',
                    'IGHD3_22','IGHD2_21', 'IGHD1_20', 'IGHD6_19',
                    'IGHD5_18', 'IGHD4_17','IGHD3_16', 'IGHD2_15', 'IGHD1_14',
                    'IGHD6_13', 'IGHD5_12','IGHD4_11', 'IGHD3_10', 'IGHD3_9',
                    'IGHD2_8', 'IGHD1_7','IGHD6_6' , 'IGHD5_5', 'IGHD4_4',
                    'IGHD3_3', 'IGHD2_2', 'IGHD1_1','IGHD7_27',
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
                    'IGLC7','IGLL5')
  
  # Need segment ranges start + stop 
  segment.ranges <- vdj.segments.list[['IGH']] %>%
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
  
  class.switch.genes <- c('IGHA2','IGHE','IGHG4','IGHG2','IGHA1',
                          'IGHG1','IGHG3','IGHM')
  segment.ranges <- segment.ranges %>%
    dplyr::mutate(segType = ifelse(segName %in% class.switch.genes,
                                   'CS',segType)) %>%
    dplyr::mutate(next_segType= ifelse(segName %in% class.switch.genes[-8],
                                       'CS',next_segType)) %>%
    dplyr::mutate(prev_segType = ifelse(segName %in% c(class.switch.genes,'IGHJ6'),
                                        'CS',next_segType))
  
  v_seg_loc <- which(segment.ranges$next_segType == 'V')
  new_seg_name <- paste0(segment.ranges$segName[v_seg_loc[1]],'_',
                         segment.ranges$segName[v_seg_loc[1] + 1])
  segment.ranges$segName[v_seg_loc[1]] <- new_seg_name
  segment.ranges <- segment.ranges[-c(v_seg_loc[1] + 1), ]
  segment.ranges <- segment.ranges %>%
    dplyr::select(segName, start = start2, end = end2) %>%
    dplyr::mutate(segName = as.character(segName))
  
  # Does segments make sense here????
  segs_to_remove <- setdiff(segment.ranges$segName, model.fit$segment)
  
  segment.ranges <- segment.ranges %>%
    filter(!segName %in% segs_to_remove) %>%
    mutate(range_check = end+1 == lead(start, default = TRUE)) %>%
    mutate(range_check = ifelse(segName == 'IGHG3', TRUE, range_check))
  range_check_loc <- which(!segment.ranges$range_check)
  range_check_loc_v <- intersect(which(grepl('IGHV',segment.ranges$segName)), range_check_loc)
  range_check_loc_j <- intersect(which(grepl('IGHJ',segment.ranges$segName)), range_check_loc)
  if(length(range_check_loc_v) > 0){
    for(i1 in 1:length(range_check_loc_v)){
      segment.ranges$end[range_check_loc_v] <- segment.ranges$start[range_check_loc_v + 1] -1
    }
  }
  if(length(range_check_loc_j) > 0){
    for(i1 in 1:length(range_check_loc_j)){
      segment.ranges$start[range_check_loc_j + 1] <- segment.ranges$end[range_check_loc_j] +1
    }
  }
  segment.ranges <- segment.ranges %>%
    dplyr::select(segName, start, end)
  
  
  # Add segments to df
  for(i in seq_len(dim(segment.ranges)[1])){
    input_cov_gc <- create_col_seg(input_cov_gc,
                                 start_pos = segment.ranges$start[i],
                                 end_pos = segment.ranges$end[i],
                                 col.name = segment.ranges$segName[i])
  }
  
  # Fit actual model to data
  input_cov_gc$logR_fit <- 0
  for(i in seq_len(length(model.fit$segment))){
    seg_input <- model.fit$segment[i]
    logR_input <- model.fit$logR[i]
    input_cov_gc$logR_fit <- input_cov_gc$logR_fit + input_cov_gc[[seg_input]] * logR_input
  }
  
  
  # Calculate corresponding reads
  input_cov_gc$reads_fit <- (2^input_cov_gc$logR_fit)*n1
  # Calculate read difference to return
  input_cov_gc$reads_diff <- input_cov_gc$reads - input_cov_gc$reads_fit
  # Also calculate fraction of the difference
  input_cov_gc$reads_diff_frac <- input_cov_gc$reads_diff/input_cov_gc$reads_fit
  input_cov_gc$reads_norm <- median(input_cov_gc$reads) + (input_cov_gc$reads_diff_frac * median(input_cov_gc$reads))
  
  input_cov_gc_out <- input_cov_gc %>%
    dplyr::select(pos,reads, Ratio, logR_fit, reads_fit, reads_diff,reads_diff_frac, reads_norm)
  return(input_cov_gc_out)
  
}


calc_haplotype_regions_de_novo_diploid <- function(input_cov, kb_len_threshold = 5){
  # Assumes diploid (run on germline sample if possible)
  cov_df_segment <- input_cov %>%
    dplyr::mutate(reads_norm = reads/median(input_cov$reads, na.rm = TRUE)) %>%
    dplyr::mutate(reads_norm = zoo::rollmedian(reads_norm, k = 1001, fill = NA)) %>%
    dplyr::filter(!is.na(reads_norm)) %>%
    dplyr::mutate(reads_seg = round(2*reads_norm)/2)
  
  cov_df_segment_summary <- cov_df_segment %>%
    dplyr::mutate(kb_bin = pos %/% 1000) %>%
    dplyr::group_by(kb_bin) %>%
    dplyr::summarise(reads_seg_average = mean(reads_seg)) %>%
    dplyr::mutate(reads_seg_average = round(2*reads_seg_average)/2) %>%
    dplyr::mutate(reads_seg_average_roll = zoo::rollmedian(reads_seg_average,
                                                           k = 11, fill = NA))
  
  ranges_summary <- cov_df_segment_summary %>%
    dplyr::filter(!is.na(reads_seg_average_roll)) %>%
    dplyr::mutate(change = ifelse(lag(reads_seg_average_roll, default = 0) != reads_seg_average_roll, 1,0)) %>%
    dplyr::mutate(change2 = cumsum(change)) %>%
    dplyr::group_by(change2) %>%
    dplyr::summarise(start_kb = dplyr::first(kb_bin),
                     end_kb = dplyr::last(kb_bin),
                     reads_seg_average_roll = dplyr::first(reads_seg_average_roll)) %>%
    dplyr::mutate(kb_length = 1 + end_kb - start_kb)
  
  ranges_summary2 <- ranges_summary %>%
    dplyr::mutate(previous_and_current_greater_1 = lag(reads_seg_average_roll > 1, default = FALSE) & reads_seg_average_roll > 1) %>%
    dplyr::mutate(next_and_current_greater_1 = lead(reads_seg_average_roll > 1, default = FALSE) & reads_seg_average_roll > 1) %>%
    dplyr::mutate(greater_1_group = previous_and_current_greater_1 | next_and_current_greater_1) %>%
    dplyr::mutate(greater_1_group_clust = cumsum(greater_1_group & !lag(greater_1_group, default = FALSE))) %>%
    dplyr::mutate(greater_1_group_clust = ifelse(greater_1_group, greater_1_group_clust, NA)) %>%
    dplyr::select(-greater_1_group, -previous_and_current_greater_1, -next_and_current_greater_1) %>%
    dplyr::mutate(previous_and_current_less_1 = lag(reads_seg_average_roll < 1, default = FALSE) & reads_seg_average_roll < 1) %>%
    dplyr::mutate(next_and_current_less_1 = lead(reads_seg_average_roll < 1, default = FALSE) & reads_seg_average_roll < 1) %>%
    dplyr::mutate(less_1_group = previous_and_current_less_1 | next_and_current_less_1) %>%
    dplyr::mutate(less_1_group_clust = cumsum(less_1_group & !lag(less_1_group, default = FALSE))) %>%
    dplyr::mutate(less_1_group_clust = ifelse(less_1_group, less_1_group_clust, NA)) %>%
    dplyr::select(-less_1_group, -previous_and_current_less_1, -next_and_current_less_1) 
    
  ranges_summary2_single_ranges <- ranges_summary2 %>%
    dplyr::filter(is.na(greater_1_group_clust) & is.na(less_1_group_clust)) %>%
    dplyr::select(start_kb, end_kb, reads_seg_average_roll, kb_length) %>%
    dplyr::mutate(kb_group_length = kb_length)
  
  tmp_df1 <- ranges_summary2 %>%
    dplyr::filter(!is.na(greater_1_group_clust)) %>%
    dplyr::group_by(greater_1_group_clust) %>%
    dplyr::summarise(kb_group_length = sum(kb_length)) %>%
    dplyr::select(greater_1_group_clust, kb_group_length)
  
  ranges_summary2_greater_clust_ranges <- ranges_summary2 %>%
    dplyr::filter(!is.na(greater_1_group_clust)) %>%
    dplyr::left_join(tmp_df1, 'greater_1_group_clust') %>%
    dplyr::select(start_kb, end_kb, reads_seg_average_roll, kb_length, kb_group_length)
  
  tmp_df2 <- ranges_summary2 %>%
    dplyr::filter(!is.na(less_1_group_clust)) %>%
    dplyr::group_by(less_1_group_clust) %>%
    dplyr::summarise(kb_group_length = sum(kb_length)) %>%
    dplyr::select(less_1_group_clust, kb_group_length)
  
  ranges_summary2_less_clust_ranges <- ranges_summary2 %>%
    dplyr::filter(!is.na(less_1_group_clust)) %>%
    dplyr::left_join(tmp_df2, 'less_1_group_clust') %>%
    dplyr::select(start_kb, end_kb, reads_seg_average_roll, kb_length, kb_group_length)
  
  ranges_summary3 <- rbind(ranges_summary2_single_ranges, 
                           ranges_summary2_greater_clust_ranges) %>%
    rbind(ranges_summary2_less_clust_ranges)
  
  # Identify ranges where kb_length > 10 + not 1
  ranges_to_normalise <- ranges_summary3 %>%
    dplyr::filter(reads_seg_average_roll != 1) %>%
    dplyr::filter(kb_group_length >= kb_len_threshold)%>%
    dplyr::rename(reads_seg = reads_seg_average_roll)
  return(ranges_to_normalise)
}

calc_haplotype_regions_de_novo_diploid_noroll <- function(input_cov, kb_len_threshold = 5){
  # Assumes diploid (run on germline sample if possible)
  cov_df_segment <- input_cov %>%
    dplyr::mutate(reads_norm = reads/median(input_cov$reads, na.rm = TRUE)) %>%
    dplyr::mutate(reads_norm = zoo::rollmedian(reads_norm, k = 1001, fill = NA)) %>%
    dplyr::filter(!is.na(reads_norm)) %>%
    dplyr::mutate(reads_seg = round(2*reads_norm)/2)
  
  cov_df_segment_summary <- cov_df_segment %>%
    dplyr::mutate(kb_bin = pos %/% 1000) %>%
    dplyr::group_by(kb_bin) %>%
    dplyr::summarise(reads_seg_average = mean(reads_seg)) %>%
    dplyr::mutate(reads_seg_average = round(2*reads_seg_average)/2)
  
  ranges_summary <- cov_df_segment_summary %>%
    dplyr::filter(!is.na(reads_seg_average)) %>%
    dplyr::mutate(change = ifelse(lag(reads_seg_average, default = 0) != reads_seg_average, 1,0)) %>%
    dplyr::mutate(change2 = cumsum(change)) %>%
    dplyr::group_by(change2) %>%
    dplyr::summarise(start_kb = dplyr::first(kb_bin),
                     end_kb = dplyr::last(kb_bin),
                     reads_seg_average = dplyr::first(reads_seg_average)) %>%
    dplyr::mutate(kb_length = 1 + end_kb - start_kb)
  
  ranges_summary2 <- ranges_summary %>%
    dplyr::mutate(previous_and_current_greater_1 = lag(reads_seg_average > 1, default = FALSE) & reads_seg_average > 1) %>%
    dplyr::mutate(next_and_current_greater_1 = lead(reads_seg_average > 1, default = FALSE) & reads_seg_average > 1) %>%
    dplyr::mutate(greater_1_group = previous_and_current_greater_1 | next_and_current_greater_1) %>%
    dplyr::mutate(greater_1_group_clust = cumsum(greater_1_group & !lag(greater_1_group, default = FALSE))) %>%
    dplyr::mutate(greater_1_group_clust = ifelse(greater_1_group, greater_1_group_clust, NA)) %>%
    dplyr::select(-greater_1_group, -previous_and_current_greater_1, -next_and_current_greater_1) %>%
    dplyr::mutate(previous_and_current_less_1 = lag(reads_seg_average < 1, default = FALSE) & reads_seg_average < 1) %>%
    dplyr::mutate(next_and_current_less_1 = lead(reads_seg_average < 1, default = FALSE) & reads_seg_average < 1) %>%
    dplyr::mutate(less_1_group = previous_and_current_less_1 | next_and_current_less_1) %>%
    dplyr::mutate(less_1_group_clust = cumsum(less_1_group & !lag(less_1_group, default = FALSE))) %>%
    dplyr::mutate(less_1_group_clust = ifelse(less_1_group, less_1_group_clust, NA)) %>%
    dplyr::select(-less_1_group, -previous_and_current_less_1, -next_and_current_less_1) 
  
  ranges_summary2_single_ranges <- ranges_summary2 %>%
    dplyr::filter(is.na(greater_1_group_clust) & is.na(less_1_group_clust)) %>%
    dplyr::select(start_kb, end_kb, reads_seg_average, kb_length) %>%
    dplyr::mutate(kb_group_length = kb_length)
  
  tmp_df1 <- ranges_summary2 %>%
    dplyr::filter(!is.na(greater_1_group_clust)) %>%
    dplyr::group_by(greater_1_group_clust) %>%
    dplyr::summarise(kb_group_length = sum(kb_length)) %>%
    dplyr::select(greater_1_group_clust, kb_group_length)
  
  ranges_summary2_greater_clust_ranges <- ranges_summary2 %>%
    dplyr::filter(!is.na(greater_1_group_clust)) %>%
    dplyr::left_join(tmp_df1, 'greater_1_group_clust') %>%
    dplyr::select(start_kb, end_kb, reads_seg_average, kb_length, kb_group_length)
  
  tmp_df2 <- ranges_summary2 %>%
    dplyr::filter(!is.na(less_1_group_clust)) %>%
    dplyr::group_by(less_1_group_clust) %>%
    dplyr::summarise(kb_group_length = sum(kb_length)) %>%
    dplyr::select(less_1_group_clust, kb_group_length)
  
  ranges_summary2_less_clust_ranges <- ranges_summary2 %>%
    dplyr::filter(!is.na(less_1_group_clust)) %>%
    dplyr::left_join(tmp_df2, 'less_1_group_clust') %>%
    dplyr::select(start_kb, end_kb, reads_seg_average, kb_length, kb_group_length)
  
  ranges_summary3 <- rbind(ranges_summary2_single_ranges, 
                           ranges_summary2_greater_clust_ranges) %>%
    rbind(ranges_summary2_less_clust_ranges)
  
  # Identify ranges where kb_length > 10 + not 1
  ranges_to_normalise <- ranges_summary3 %>%
    dplyr::filter(reads_seg_average != 1) %>%
    dplyr::filter(kb_length >= kb_len_threshold) %>%
    dplyr::rename(reads_seg = reads_seg_average)
  return(ranges_to_normalise)
}


calc_haplotype_regions_de_novo_highBcell_nondiploid <- function(input_cov, kb_len_threshold = 5,
                                                                round_solution = TRUE){
  V_start <- (105860500 %/% 1000) + 1
  CS_end <- 105860500 %/% 1000

  # Assumes diploid (run on germline sample if possible)
  cov_df_segment <- input_cov %>%
    dplyr::mutate(reads_norm = reads/median(tail(input_cov$reads,n = 20000), na.rm = TRUE)) %>%
    dplyr::mutate(reads_norm = zoo::rollmedian(reads_norm, k = 1001, fill = NA)) %>%
    dplyr::filter(!is.na(reads_norm))
    
  # Calculate 1kb bins coverage value
  cov_df_segment_summary <- cov_df_segment %>%
    dplyr::mutate(kb_bin = pos %/% 1000) %>%
    dplyr::group_by(kb_bin) %>%
    dplyr::summarise(reads_seg_average = mean(reads_norm)) 
  
  cov_df_segment_summary_v_loc <- cov_df_segment_summary %>%
    filter(kb_bin > V_start)
  cov_df_segment_summary_cs_loc <- cov_df_segment_summary %>%
    filter(kb_bin < CS_end) 
  
  # Classify potential duplication or deletion events:
  cov_df_segment_summary_v_loc$deletion_flag <- 0
  cov_df_segment_summary_v_loc$duplication_flag <- 0
  for(i in 10:(dim(cov_df_segment_summary_v_loc)[1] - 10)){
    rev.idx <- rev(1:(dim(cov_df_segment_summary_v_loc)[1]))
    
    # Get maximum of 100 reads
    
    forward.reads <- cov_df_segment_summary_v_loc[seq(i),] %>%
      filter(deletion_flag == 0) %>%
      filter(duplication_flag == 0) %>%
      select(reads_seg_average) %>%
      `[[`(1) 
    if(length(forward.reads) > 50){
      forward.reads <- tail(forward.reads, 50)
    }
    
    forward.median <- median(forward.reads)
    forward.sd <- sd(forward.reads)
    below_thresh <- forward.median - 2*forward.sd 
    below_thresh_loc <- which(cov_df_segment_summary_v_loc$reads_seg_average < below_thresh)
    below_thresh_loc <- setdiff(below_thresh_loc, seq(i))
    cov_df_segment_summary_v_loc$deletion_flag[below_thresh_loc] <- 1
    
    backward.reads  <-  cov_df_segment_summary_v_loc[rev.idx[seq(i)],] %>%
      filter(deletion_flag == 0) %>%
      filter(duplication_flag == 0) %>%
      select(reads_seg_average) %>%
      `[[`(1)
    if(length(backward.reads) > 50){
      backward.reads <- tail(backward.reads, 50)
    }
    backwards.median <- median(backward.reads)
    backwards.sd <- sd(backward.reads)
    above_thresh <- backwards.median + 2*backwards.sd 
    above_thresh_loc <- which(cov_df_segment_summary_v_loc$reads_seg_average > above_thresh)
    above_thresh_loc <- setdiff(above_thresh_loc, rev.idx[seq(i)])
    cov_df_segment_summary_v_loc$duplication_flag[above_thresh_loc] <- 1
  }
  
  cov_df_segment_summary_cs_loc$deletion_flag <- 0
  cov_df_segment_summary_cs_loc$duplication_flag <- 0
  
  for(i in 10:(dim(cov_df_segment_summary_cs_loc)[1] - 10)){
    rev.idx <- rev(1:(dim(cov_df_segment_summary_cs_loc)[1]))
    
    # Get maximum of 100 reads
    
    forward.reads <- cov_df_segment_summary_cs_loc[seq(i),] %>%
      filter(deletion_flag == 0) %>%
      filter(duplication_flag == 0) %>%
      select(reads_seg_average) %>%
      `[[`(1) 
    if(length(forward.reads) > 25){
      forward.reads <- tail(forward.reads, 25)
    }
    
    forward.median <- median(forward.reads)
    forward.sd <- sd(forward.reads)
    above_thresh <- forward.median + 2*forward.sd 
    above_thresh_loc <- which(cov_df_segment_summary_cs_loc$reads_seg_average > above_thresh)
    above_thresh_loc <- setdiff(above_thresh_loc, seq(i))
    cov_df_segment_summary_cs_loc$duplication_flag[above_thresh_loc] <- 1
    
    backward.reads  <-  cov_df_segment_summary_cs_loc[rev.idx[seq(i)],] %>%
      filter(deletion_flag == 0) %>%
      filter(duplication_flag == 0) %>%
      select(reads_seg_average) %>%
      `[[`(1)
    if(length(backward.reads) > 25){
      backward.reads <- tail(backward.reads, 25)
    }
    backwards.median <- median(backward.reads)
    backwards.sd <- sd(backward.reads)
    below_thresh <- backwards.median - 2*backwards.sd 
    below_thresh_loc <- which(cov_df_segment_summary_cs_loc$reads_seg_average < below_thresh)
    below_thresh_loc <- setdiff(below_thresh_loc, rev.idx[seq(i)])
    cov_df_segment_summary_cs_loc$deletion_flag[below_thresh_loc] <- 1
  }
  
  ranges_summary_vloc <- cov_df_segment_summary_v_loc  %>%
    mutate(change = deletion_flag | duplication_flag) %>%
    mutate(change_state = change != lag(change, default = FALSE) |
             deletion_flag != lag(deletion_flag, default = 0) |
             duplication_flag != lag(duplication_flag, default = 0)) %>%
    mutate(change2 = cumsum(change_state)) %>%
    dplyr::group_by(change2) %>%
    dplyr::summarise(start_kb = dplyr::first(kb_bin),
                     end_kb = dplyr::last(kb_bin),
                     reads_seg_average_mean = mean(reads_seg_average),
                     change = unique(ifelse(deletion_flag == 1, 'deletion',
                                            ifelse(duplication_flag == 1, 'duplication','no_change')))) %>%
    dplyr::mutate(kb_length = 1 + end_kb - start_kb) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(reads_seg_average_baseline = ifelse(change != 'no_change', 
                                                      (lag(reads_seg_average_mean) + lead(reads_seg_average_mean))/2, NA)) %>%
    dplyr::mutate(reads_seg_average_roll = reads_seg_average_mean/reads_seg_average_baseline) %>%
    dplyr::select(start_kb, end_kb, reads_seg_average_roll, kb_length) 
  
  ranges_summary_csloc <- cov_df_segment_summary_cs_loc  %>%
    mutate(change = deletion_flag | duplication_flag) %>%
    mutate(change_state = change != lag(change, default = FALSE) |
             deletion_flag != lag(deletion_flag, default = 0) |
             duplication_flag != lag(duplication_flag, default = 0)) %>%
    mutate(change2 = cumsum(change_state)) %>%
    dplyr::group_by(change2) %>%
    dplyr::summarise(start_kb = dplyr::first(kb_bin),
                     end_kb = dplyr::last(kb_bin),
                     reads_seg_average_mean = mean(reads_seg_average),
                     change = unique(ifelse(deletion_flag == 1, 'deletion',
                                            ifelse(duplication_flag == 1, 'duplication','no_change')))) %>%
    dplyr::mutate(kb_length = 1 + end_kb - start_kb) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(reads_seg_average_baseline = ifelse(change != 'no_change', 
                                                      (lag(reads_seg_average_mean) + lead(reads_seg_average_mean))/2, NA)) %>%
    dplyr::mutate(reads_seg_average_roll = reads_seg_average_mean/reads_seg_average_baseline) %>%
    dplyr::select(start_kb, end_kb, reads_seg_average_roll, kb_length) 
  
  if(round_solution){
    ranges_to_normalise <- ranges_summary_csloc %>%
      rbind(ranges_summary_vloc) %>%
      dplyr::mutate(reads_seg_average_roll = round(2*reads_seg_average_roll)/2) %>%
      dplyr::mutate(reads_seg_average_roll = ifelse(is.na(reads_seg_average_roll), 1, reads_seg_average_roll)) %>%
      dplyr::filter(reads_seg_average_roll != 1) %>%
      dplyr::filter(kb_length >= kb_len_threshold)%>%
      dplyr::rename(reads_seg = reads_seg_average_roll)
  }else{
    ranges_to_normalise <- ranges_summary_csloc %>%
      rbind(ranges_summary_vloc) %>%
      dplyr::mutate(reads_seg_average_roll = ifelse(is.na(reads_seg_average_roll), 1, reads_seg_average_roll)) %>%
      dplyr::filter(reads_seg_average_roll != 1) %>%
      dplyr::filter(kb_length >= kb_len_threshold) %>%
      dplyr::rename(reads_seg = reads_seg_average_roll)
    
  }

  
  

  
  
  
  return(ranges_to_normalise)
}

IGH_haplotype_norm_for_region_diploid <- function(input_cov, region_df){
  # Region has already been defined
  ranges_to_normalise_long <- region_df %>%
    dplyr::rowwise() %>%
    dplyr::transmute(
      pos = list(seq(start_kb * 1000, end_kb * 1000 + 999)),
      value = reads_seg
    ) %>%
    tidyr::unnest(cols = c(pos))
  
  # If kb_bin in ranges to normalise update based on reads_seg_average_roll value
  
  update_cov <- input_cov %>% 
    dplyr::mutate(kb_bin = pos %/% 1000) %>%
    dplyr::left_join(ranges_to_normalise_long, 'pos') %>%
    dplyr::mutate(value = ifelse(is.na(value),1, value)) %>%
    dplyr::filter(value != 0) %>%
    dplyr::mutate(update_reads = reads*(1/value)) %>%
    dplyr::select(pos, reads = update_reads)
  
  return(update_cov)
  
}

IGH_haplotype_norm_for_region_tumour <- function(input_cov, region_df){
  # We do not assume any knowledge on the tumour purity and allele specific copy-number but for each 
  # region will calculate the median change in coverage and use that
  # as the basis of the normalisation values. 
  norm_values <- sapply(seq_len(dim(region_df)[1]),
                        function(x) {
                          calc_norm_for_region(input_cov, region_df$start_kb[x],
                                               region_df$end_kb[x])})
  region_df$norm <- norm_values 
  # Put limits to prevent extreme cases
  region_df$norm <- ifelse(region_df$norm  > 3, 3,
                           ifelse(region_df$norm < 1/3, 1/3, region_df$norm))
  
  # Norm and values must be consistent
  # e.g. if norm > 1, value must be less than 1
  region_df <- region_df %>%
    dplyr::mutate(norm = ifelse(norm > 1 & reads_seg > 1, 1, norm )) %>%
    dplyr::mutate(norm = ifelse(norm < 1 & reads_seg < 1, 1, norm))
  
  # Region has already been defined
  ranges_to_normalise_long <- region_df %>%
    dplyr::rowwise() %>%
    dplyr::transmute(
      pos = list(seq(start_kb * 1000, end_kb * 1000 + 999)),
      value = reads_seg,
      norm = norm
    ) %>%
    tidyr::unnest(cols = c(pos))
  
  
  # If kb_bin in ranges to normalise update based on reads_seg_average_roll value
  
  update_cov <- input_cov %>% 
    dplyr::mutate(kb_bin = pos %/% 1000) %>%
    dplyr::left_join(ranges_to_normalise_long, 'pos') %>%
    dplyr::mutate(value = ifelse(is.na(value), 1,  value)) %>%
    dplyr::filter(value != 0) %>%
    dplyr::mutate(update_reads = ifelse(value == 1, reads, reads*norm)) %>%
    dplyr::select(pos, reads = update_reads)
  
  return(update_cov)
  
}

calc_norm_for_region <- function(input_cov, start_pos, end_pos, kb_range = 2500){
  # Median coverage in region:
  baseline_median <- input_cov %>%
    filter(pos < 105588395 | pos > 106779812)  %>%
    select(reads) %>% `[[`(1) %>% median(na.rm = TRUE)
  
  start_pos <- start_pos*1000
  end_pos <- end_pos*1000 + 999
  region_median <- input_cov %>%
    filter(pos >= start_pos) %>%
    filter(pos <= end_pos) %>%
    select(reads) %>% `[[`(1) %>% median(na.rm = TRUE)
  # Median coverage in kb_range before region:
  # pre_region_median <- input_cov %>%
  #   filter(pos >= start_pos - kb_range) %>%
  #   filter(pos <= start_pos) %>%
  #   select(reads) %>% `[[`(1) %>% median(na.rm = TRUE)
  # 
  # # Median coverage in kb_range after region:
  # post_region_median <- input_cov %>%
  #   filter(pos >= end_pos + 999) %>%
  #   filter(pos <= end_pos + 1000 + kb_range) %>%
  #   select(reads) %>% `[[`(1) %>% median(na.rm = TRUE)
  
  # Before vs region:
  region_change <- baseline_median/region_median
  # before_vs_region <- pre_region_median/region_median
  # # After vs regions:
  # post_vs_region <- post_region_median/region_median
  #  
  # region_change = (before_vs_region + post_vs_region)/2
  # if(is.na(before_vs_region) & !is.na(post_vs_region)){
  #   region_change = post_vs_region
  # }
  # if(is.na(post_vs_region) & !is.na(before_vs_region)){
  #   region_change = before_vs_region
  # }
  # if(is.na(post_vs_region) & is.na(before_vs_region)){
  #   region_change = 1
  # }
  norm.value = region_change
  
  return(norm.value)
  
}