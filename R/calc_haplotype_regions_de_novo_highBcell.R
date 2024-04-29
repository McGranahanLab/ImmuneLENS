#' Calculate predicted IGH germline haplotype assuming possible high B cell content
#'
#' @param input_cov data frame of positions and coverage values for IGH loci
#' @param kb_len_threshold threshold for length of genomic regions used in normalisation method (default = 5kb)
#' @param round_solution default = TRUE
#' @return  regions used for later normalisation
#' @name calc_haplotype_regions_de_novo_highBcell

calc_haplotype_regions_de_novo_highBcell <- function(input_cov, kb_len_threshold = 5,
                                                     round_solution = TRUE){
  reads <- pos <- reads_norm <- reads_seg <- kb_bin <- reads_seg_average <- reads_seg_average_roll <- change <- change2 <- end_kb <- start_kb <- NULL
  deletion_flag <- duplication_flag <- change_state <-  reads_seg_average_mean <- reads_seg_average_roll <- reads_seg_average_baseline <- kb_length <- germline_reads <-  NULL
  
  V_start <- (105860500 %/% 1000) + 1
  CS_end <- 105860500 %/% 1000
  
  # Assumes diploid (run on germline sample if possible)
  cov_df_segment <- input_cov %>%
    dplyr::mutate(reads_norm = reads/median(utils::tail(input_cov$reads,n = 20000), na.rm = TRUE)) %>%
    dplyr::mutate(reads_norm = zoo::rollmedian(reads_norm, k = 1001, fill = NA)) %>%
    dplyr::filter(!is.na(reads_norm))
  
  # Calculate 1kb bins coverage value
  cov_df_segment_summary <- cov_df_segment %>%
    dplyr::mutate(kb_bin = pos %/% 1000) %>%
    dplyr::group_by(kb_bin) %>%
    dplyr::summarise(reads_seg_average = mean(reads_norm)) 
  
  cov_df_segment_summary_v_loc <- cov_df_segment_summary %>%
    dplyr::filter(kb_bin > V_start)
  cov_df_segment_summary_cs_loc <- cov_df_segment_summary %>%
    dplyr::filter(kb_bin < CS_end) 
  
  # Classify potential duplication or deletion events:
  cov_df_segment_summary_v_loc$deletion_flag <- 0
  cov_df_segment_summary_v_loc$duplication_flag <- 0
  for(i in 10:(dim(cov_df_segment_summary_v_loc)[1] - 10)){
    rev.idx <- rev(1:(dim(cov_df_segment_summary_v_loc)[1]))
    
    # Get maximum of 100 reads
    
    forward.reads <- cov_df_segment_summary_v_loc[seq(i),] %>%
      dplyr::filter(deletion_flag == 0) %>%
      dplyr::filter(duplication_flag == 0) %>%
      dplyr::select(reads_seg_average) %>%
      `[[`(1) 
    if(length(forward.reads) > 50){
      forward.reads <- utils::tail(forward.reads, 50)
    }
    
    forward.median <- median(forward.reads)
    forward.sd <- sd(forward.reads)
    below_thresh <- forward.median - 2*forward.sd 
    below_thresh_loc <- which(cov_df_segment_summary_v_loc$reads_seg_average < below_thresh)
    below_thresh_loc <- setdiff(below_thresh_loc, seq(i))
    cov_df_segment_summary_v_loc$deletion_flag[below_thresh_loc] <- 1
    
    backward.reads  <-  cov_df_segment_summary_v_loc[rev.idx[seq(i)],] %>%
      dplyr::filter(deletion_flag == 0) %>%
      dplyr::filter(duplication_flag == 0) %>%
      dplyr::select(reads_seg_average) %>%
      `[[`(1)
    if(length(backward.reads) > 50){
      backward.reads <-  utils::tail(backward.reads, 50)
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
      dplyr::filter(deletion_flag == 0) %>%
      dplyr::filter(duplication_flag == 0) %>%
      dplyr::select(reads_seg_average) %>%
      `[[`(1) 
    if(length(forward.reads) > 25){
      forward.reads <-  utils::tail(forward.reads, 25)
    }
    
    forward.median <- median(forward.reads)
    forward.sd <- sd(forward.reads)
    above_thresh <- forward.median + 2*forward.sd 
    above_thresh_loc <- which(cov_df_segment_summary_cs_loc$reads_seg_average > above_thresh)
    above_thresh_loc <- setdiff(above_thresh_loc, seq(i))
    cov_df_segment_summary_cs_loc$duplication_flag[above_thresh_loc] <- 1
    
    backward.reads  <-  cov_df_segment_summary_cs_loc[rev.idx[seq(i)],] %>%
      dplyr::filter(deletion_flag == 0) %>%
      dplyr::filter(duplication_flag == 0) %>%
      dplyr::select(reads_seg_average) %>%
      `[[`(1)
    if(length(backward.reads) > 25){
      backward.reads <-  utils::tail(backward.reads, 25)
    }
    backwards.median <- median(backward.reads)
    backwards.sd <- sd(backward.reads)
    below_thresh <- backwards.median - 2*backwards.sd 
    below_thresh_loc <- which(cov_df_segment_summary_cs_loc$reads_seg_average < below_thresh)
    below_thresh_loc <- setdiff(below_thresh_loc, rev.idx[seq(i)])
    cov_df_segment_summary_cs_loc$deletion_flag[below_thresh_loc] <- 1
  }
  
  ranges_summary_vloc <- cov_df_segment_summary_v_loc  %>%
    dplyr::mutate(change = deletion_flag | duplication_flag) %>%
    dplyr::mutate(change_state = change != dplyr::lag(change, default = FALSE) |
             deletion_flag != dplyr::lag(deletion_flag, default = 0) |
             duplication_flag != dplyr::lag(duplication_flag, default = 0)) %>%
    dplyr::mutate(change2 = cumsum(change_state)) %>%
    dplyr::group_by(change2) %>%
    dplyr::summarise(start_kb = dplyr::first(kb_bin),
                     end_kb = dplyr::last(kb_bin),
                     reads_seg_average_mean = mean(reads_seg_average),
                     change = unique(ifelse(deletion_flag == 1, 'deletion',
                                            ifelse(duplication_flag == 1, 'duplication','no_change')))) %>%
    dplyr::mutate(kb_length = 1 + end_kb - start_kb) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(reads_seg_average_baseline = ifelse(change != 'no_change', 
                                                      (dplyr::lag(reads_seg_average_mean) + dplyr::lead(reads_seg_average_mean))/2, NA)) %>%
    dplyr::mutate(reads_seg_average_roll = reads_seg_average_mean/reads_seg_average_baseline) %>%
    dplyr::select(start_kb, end_kb, reads_seg_average_roll, kb_length) 
  
  ranges_summary_csloc <- cov_df_segment_summary_cs_loc  %>%
    dplyr::mutate(change = deletion_flag | duplication_flag) %>%
    dplyr::mutate(change_state = change != dplyr::lag(change, default = FALSE) |
             deletion_flag != dplyr::lag(deletion_flag, default = 0) |
             duplication_flag != dplyr::lag(duplication_flag, default = 0)) %>%
    dplyr::mutate(change2 = cumsum(change_state)) %>%
    dplyr::group_by(change2) %>%
    dplyr::summarise(start_kb = dplyr::first(kb_bin),
                     end_kb = dplyr::last(kb_bin),
                     reads_seg_average_mean = mean(reads_seg_average),
                     change = unique(ifelse(deletion_flag == 1, 'deletion',
                                            ifelse(duplication_flag == 1, 'duplication','no_change')))) %>%
    dplyr::mutate(kb_length = 1 + end_kb - start_kb) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(reads_seg_average_baseline = ifelse(change != 'no_change', 
                                                      (dplyr::lag(reads_seg_average_mean) + dplyr::lead(reads_seg_average_mean))/2, NA)) %>%
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
