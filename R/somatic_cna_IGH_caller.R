#' Calculate predicted regions containing somatic copy number alterations in IGH locus
#'
#' @param germline_cov data frame of positions and coverage values for IGH loci of germline sample
#' @param tumour_cov data frame of positions and coverage values for IGH loci of tumour sample
#' @param gc_correct logical, should the normalisation be GC corrected? (default = TRUE)
#' @param run_breakpoint_check logical, should the normalisation be checked for breakpoints? (default = TRUE)
#' @param rpart_cp complexity parameter for rpart model (default = 0.01)
#' @return Data frame containing regions of predicted somatic copy number alterations
#' @name somatic_cna_IGH_caller

somatic_cna_IGH_caller <- function(germline_cov, tumour_cov, gc_correct = TRUE, run_breakpoint_check = FALSE, 
                                   rpart_cp = 0.01){
  reads <- germline_reads <- tumour_reads <- pos <- kb_bin <- logR <- start <- end <- start_kb <- end_kb <- check_overlap_both_sides <- ratio <- NULL
  
  # Blacklist most focal part of IGH cs and VDJ region
  # CS focal: 105712500 - 105860500 (IGHA1 - IGHM) 105712-105860
  # VDJ focal: 105865458 - 105939756 (J - V) 105865-105939
  germline_cov_orig <- germline_cov
  tumour_cov_orig <- tumour_cov
  if(gc_correct){
    germline_cov  <- IGH_gc_correct(germline_cov)
    tumour_cov  <- IGH_gc_correct(tumour_cov)
  } 
  
  logR_df <- germline_cov %>%
    dplyr::rename(germline_reads = reads) %>%
    dplyr::left_join(tumour_cov,'pos') %>%
    dplyr::rename(tumour_reads = reads) %>%
    dplyr::mutate(germline_reads = germline_reads/median(germline_cov$reads, na.rm = TRUE)) %>%
    dplyr::mutate(tumour_reads = tumour_reads/median(tumour_cov$reads, na.rm = TRUE)) %>%
    dplyr::mutate(logR = log2((tumour_reads/germline_reads) + 1)) %>%
    dplyr::mutate(kb_bin = pos %/% 1000) %>% 
    dplyr::filter(! (kb_bin >= 105712 & kb_bin <= 105860)) %>%
    dplyr::filter(! (kb_bin >= 105865 & kb_bin <= 105939)) %>%
    dplyr::group_by(kb_bin) %>%
    dplyr::summarise(germline_reads = median(germline_reads, na.rm = TRUE),
              tumour_reads = median(tumour_reads, na.rm = TRUE),
              logR = median(logR, na.rm = TRUE))
  
  logR_df <- logR_df %>%
    dplyr::filter(!is.na(logR)) %>%
    dplyr::filter(!is.infinite(logR))
  
  # Segment this into line segments - how?
  # logR_roll <- zoo::rollmedian(logR_df$logR,k = 51)
  
  tree <- rpart::rpart(logR ~ kb_bin, data=logR_df,
                       control = rpart::rpart.control(cp = rpart_cp))
  
  tree.splits.loc <- as.numeric(tree$splits[,'index'])
  tree.splits.loc <- tree.splits.loc[order(tree.splits.loc)]
  IGH.start <- vdj_seg_list[['IGH_hg38']][1,2] %/% 1000
  IGH.end <- vdj_seg_list[['IGH_hg38']][1,3] %/% 1000
  
  split.start.vec <- numeric(length(tree.splits.loc) + 1)
  split.end.vec <- numeric(length(tree.splits.loc) + 1)
  split.start.vec[1] <- IGH.start
  split.end.vec[length(tree.splits.loc) + 1] <- IGH.end
  for(i in seq_len(length(tree.splits.loc))){
    split.start.vec[i+1] <- round(tree.splits.loc)[i] 
    split.end.vec[i] <- round(tree.splits.loc)[i] - 1
  }
  sm_segment_ranges_df <- data.frame(start = split.start.vec,
                                     end = split.end.vec)
  
  
  
  # For each segment get median logR germline reads and tumour reads
  logR_summary_list <- list()
  for(i in seq_len(dim(sm_segment_ranges_df)[1])){
    logR_summary_list[[i]] <- logR_df %>%
      dplyr::filter(kb_bin >= sm_segment_ranges_df$start[i]) %>%
      dplyr::filter(kb_bin <= sm_segment_ranges_df$end[i]) %>%
      dplyr::summarise(germline_reads = median(germline_reads,na.rm = TRUE),
                tumour_reads = median(tumour_reads, na.rm = TRUE),
                logR = median(logR, na.rm = TRUE)) %>%
      dplyr::mutate(start = sm_segment_ranges_df$start[i],
             end = sm_segment_ranges_df$end[i]) %>%
      dplyr::select(start, end, germline_reads, tumour_reads, logR)
  }
  logR_summary_df <- data.table::rbindlist(logR_summary_list) %>%
    dplyr::mutate(kb_len = end - start)
  
  
  # Identify which ranges are the baseline and should not change
  baseline.loc <- which.min(abs(1 - logR_summary_df$logR))
  baseline.median <- tumour_cov_orig %>%
    dplyr::filter(pos > logR_summary_df[baseline.loc, 'start'][[1]] * 1000) %>%
    dplyr::filter(pos < logR_summary_df[baseline.loc, 'end'][[1]] * 1000) %>%
    dplyr::filter(! (pos >= 105712*1000 & pos <= 105860*1000)) %>%
    dplyr::filter(! (pos  >= 105865*1000 & pos <= 105939*1000)) %>%
    dplyr::select(reads) %>% `[[`(1) %>% median(na.rm = TRUE)
  
  # Make regions_sm_df object for tumour reads
  sm_segment_ranges_df2 <- list()
  for(i in seq_len(dim(sm_segment_ranges_df)[1])){
    test.median <- tumour_cov_orig %>%
      dplyr::filter(pos > logR_summary_df[i, 'start'][[1]] * 1000) %>%
      dplyr::filter(pos < logR_summary_df[i, 'end'][[1]] * 1000) %>%
      dplyr::filter(! (pos >= 105712*1000 & pos <= 105860*1000)) %>%
      dplyr::filter(! (pos  >= 105865*1000 & pos <= 105939*1000)) %>%
      dplyr::select(reads) %>% `[[`(1) %>% median(na.rm = TRUE)
    reads_seg <- test.median/baseline.median
    sm_segment_ranges_df2[[i]] <- data.frame(start_kb = sm_segment_ranges_df$start[i],
                                             end_kb = sm_segment_ranges_df$end[i],
                                             reads_seg = test.median/baseline.median,
                                             median_reads = test.median,
                                             kb_length =   logR_summary_df$kb_len[i],
                                             germline_reads = logR_summary_df$germline_reads[i],
                                             tumour_reads = logR_summary_df$tumour_reads[i],
                                             logR = logR_summary_df$logR[i])
  }
  sm_segment_ranges_df2 <- data.table::rbindlist(sm_segment_ranges_df2 )
  sm_segment_ranges_df2$bin_size <- 1000
  
  sm_segment_ranges_df2 <- sm_segment_ranges_df2 %>%
    dplyr::mutate(reads_seg = ifelse(abs(1-logR) < 0.05, 1, reads_seg))
  
  
  break_point_check <- sm_segment_ranges_df2 %>%
    dplyr::filter(start_kb < 105939 & end_kb > 105712) %>%
    dplyr::select(reads_seg) %>% `[[`(1) %>% `!=`(1) %>% any()
  
  if(break_point_check & run_breakpoint_check){
    segs_to_check <- sm_segment_ranges_df2 %>%
      dplyr::filter(start_kb < 105939 & end_kb > 105712)
    segs_to_check <- segs_to_check %>%
      dplyr::filter(reads_seg != 1)
    # Finally make sure that start and end are not both outside the focal regions
    segs_to_check <- segs_to_check %>% 
      dplyr::mutate(check_overlap_both_sides = start_kb < 105712 & end_kb > 105939)
    segs_to_check <- segs_to_check %>% 
      dplyr::filter(!check_overlap_both_sides)
    
    # If deletion test if it goes up suddenly within focal region
    for(i in seq_len(dim(segs_to_check)[1])){
      tmp_start <- segs_to_check$start_kb[i]
      tmp_end <- segs_to_check$end_kb[i]
      tmp_seg <- segs_to_check$reads_seg[i]
      
      if(tmp_seg == 1) { next}
      # Median of reads before the focal region
      if(tmp_start < 105712){
        tmp_median <- tumour_cov_orig %>%
          dplyr::filter(pos >= tmp_start * 1000) %>%
          dplyr::filter(pos < 105712 * 1000) %>%
          dplyr::select(reads) %>% `[[`(1) %>% median(na.rm = TRUE)
      }else{
        if(tmp_start > 105712 & tmp_end > 105939){
          tmp_median <- tumour_cov_orig %>%
            dplyr::filter(pos > 105939 * 1000) %>%
            dplyr::filter(pos <= tmp_end * 1000) %>%
            dplyr::select(reads) %>% `[[`(1) %>% median(na.rm = TRUE)
          
        }else{
          next
        }
      }
      
      tmp_start_max <- max(tmp_start, 105712)
      tmp_end_min <- min(tmp_end, 105939)
      
      tmp_ratio <- tumour_cov_orig %>%
        dplyr::filter(pos >=  tmp_start_max * 1000) %>%
        dplyr::filter(pos < tmp_end_min*1000) %>%
        dplyr::mutate(ratio = reads/tmp_median) %>%
        dplyr::mutate(kb_bin = pos %/% 1000) %>%
        dplyr::group_by(kb_bin) %>%
        dplyr::summarise(reads_average = median(ratio, na.rm = TRUE))
      
      tree2 <- rpart::rpart(reads_average~ kb_bin, data=tmp_ratio)
      tree.splits.loc2 <- as.numeric(tree2$splits[,'index'])
      tree.splits.loc2 <- tree.splits.loc2[order(tree.splits.loc2)]
      
      split.start.vec2 <- numeric(length(tree.splits.loc2) + 1)
      split.end.vec2 <- numeric(length(tree.splits.loc2) + 1)
      split.start.vec2[1] <- tmp_start_max
      split.end.vec2[length(tree.splits.loc2) + 1] <- tmp_end
      for(i1 in seq_len(length(tree.splits.loc2))){
        split.start.vec2[i1+1] <- round(tree.splits.loc2)[i1]
        split.end.vec2[i1] <- round(tree.splits.loc2)[i1] - 1
      }
      focal_segment_ranges_df <- data.frame(start = split.start.vec2,
                                            end = split.end.vec2)
      
      focal_segment_list <- list()
      for(j in seq_len(dim(focal_segment_ranges_df)[1])){
        focal_segment_list[[j]] <-  tumour_cov_orig %>%
          dplyr::filter(pos >= focal_segment_ranges_df$start[j] *1000) %>%
          dplyr::filter(pos <= focal_segment_ranges_df$end[j]*1000) %>%
          dplyr::summarise(tumour_reads = median(reads, na.rm = TRUE)) %>%
          dplyr::mutate(start = focal_segment_ranges_df$start[j],
                 end = focal_segment_ranges_df$end[j]) %>%
          dplyr::mutate(ratio = tumour_reads/tmp_median) %>%
          dplyr::select(start, end,tumour_reads, ratio)
      }
      norm_value <- 1/segs_to_check$reads_seg[i]
      focal_summary_df <- data.table::rbindlist(focal_segment_list) %>%
        dplyr::mutate(kb_len = end - start) %>%
        dplyr::mutate(tumour_reads_compared_to_norm = norm_value/ratio)
      
      # IF the norm value looks like it matches a change set this as the breakpoint!
      breakpoint_change <- abs(focal_summary_df$tumour_reads_compared_to_norm -1)
      breakpoint_change_loc <- which(breakpoint_change < 0.2)
      if(length(breakpoint_change_loc) > 0){
        # There is a breakpoint detected in the focal region
        # Check if the segment is before or after
        if(tmp_start < 105712){
          new_end <- focal_summary_df$end[breakpoint_change_loc[1]]
          # Update the somatic regions selected
          sm_segment_ranges_df2_update <- sm_segment_ranges_df2
          update_loc <- which(sm_segment_ranges_df2_update$start_kb == tmp_start)
          sm_segment_ranges_df2_update$end_kb[update_loc] <- new_end
          sm_segment_ranges_df2_update$kb_length[update_loc] <- sm_segment_ranges_df2_update$end_kb[update_loc] - sm_segment_ranges_df2_update$start_kb[update_loc]
          
          # Add in new segment (only in focal):
          new_segment <- data.frame(start_kb = new_end +1, end_kb = tmp_end, reads_seg = 1, median_reads = NA, 
                                    kb_length = tmp_end - (new_end + 1), germline_reads = NA, tumour_reads = NA,
                                    logR = NA, bin_size = 1000)
          
          sm_segment_ranges_df2_update <- rbind(sm_segment_ranges_df2_update, new_segment) %>%
            dplyr::arrange(start_kb)
          sm_segment_ranges_df2 <- sm_segment_ranges_df2_update
        }
        if(tmp_start > 105712 & tmp_end > 105939){
          new_start <- focal_summary_df$start[rev(breakpoint_change_loc)[1]]
          sm_segment_ranges_df2_update <- sm_segment_ranges_df2
          update_loc <- which(sm_segment_ranges_df2_update$start_kb == tmp_start)
          sm_segment_ranges_df2_update$start_kb[update_loc] <- new_start
          sm_segment_ranges_df2_update$kb_length[update_loc] <- sm_segment_ranges_df2_update$end_kb[update_loc] - sm_segment_ranges_df2_update$start_kb[update_loc]
          new_segment <- data.frame(start_kb = tmp_start, end_kb = new_start -1, reads_seg = 1, median_reads = NA, 
                                    kb_length = (new_start - 1) - tmp_start, germline_reads = NA, tumour_reads = NA,
                                    logR = NA, bin_size = 1000)
          sm_segment_ranges_df2_update <- rbind(sm_segment_ranges_df2_update, new_segment) %>%
            dplyr::arrange(start_kb)
          sm_segment_ranges_df2 <- sm_segment_ranges_df2_update
          
        }
      }
      
      
    }
    
  }
  
  return(sm_segment_ranges_df2)
}


