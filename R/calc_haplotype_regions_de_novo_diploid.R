#' Calculate predicted IGH germline haplotype assuming low B cell content + diploid sample
#'
#' @param input_cov data frame of positions and coverage values for IGH loci
#' @param kb_len_threshold threshold for length of genomic regions used in normalisation method (default = 5kb)
#' @param bin_size default = 1000
#' @return  regions used for later normalisation
#' @name calc_haplotype_regions_de_novo_diploid


calc_haplotype_regions_de_novo_diploid <- function(input_cov, kb_len_threshold = 5, bin_size = 1000){
  # Assumes diploid (run on germline sample if possible)
  norm.median <- median(input_cov$reads, na.rm = TRUE)
  
  for(i in 1:5){
    cov_df_segment <- input_cov %>%
      dplyr::mutate(reads_norm = reads/norm.median) %>%
      dplyr::mutate(reads_norm = zoo::rollmedian(reads_norm, k = bin_size + 1, fill = NA)) %>%
      dplyr::filter(!is.na(reads_norm)) %>%
      dplyr::mutate(reads_seg = round(2*reads_norm)/2)
    
    cov_df_segment_summary <- cov_df_segment %>%
      dplyr::mutate(kb_bin = pos %/% bin_size) %>%
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
      dplyr::rename(reads_seg = reads_seg_average_roll) %>%
      dplyr::mutate(bin_size = bin_size)
    
    # Get median coverage values based on this:
    if(dim(ranges_to_normalise)[1] > 0){ 
      ranges_to_normalise_long <- ranges_to_normalise %>%
        dplyr::rowwise() %>%
        dplyr::mutate(kb_bin = list(start_kb:end_kb)) %>%
        tidyr::unnest(kb_bin) %>%
        select(-start_kb, -end_kb) %>%
        select(kb_bin, reads_seg)  
      
      norm.median <- input_cov %>%
        dplyr::mutate(kb_bin = pos %/% bin_size) %>%
        dplyr::left_join(ranges_to_normalise_long, 'kb_bin') %>%
        dplyr::mutate(reads_seg = ifelse(is.na(reads_seg), 1, reads_seg)) %>%
        dplyr::filter(reads_seg == 1) %>%
        dplyr::select(reads) %>% `[[`(1) %>% median(na.rm = TRUE)
      
    }else{
      norm.median <- input_cov %>% 
        dplyr::mutate(kb_bin = pos %/% bin_size) %>%
        dplyr::mutate(reads_seg = 1) %>%
        dplyr::filter(reads_seg == 1) %>%
        dplyr::select(reads) %>% `[[`(1) %>% median(na.rm = TRUE)
    }
  }
  
  return(ranges_to_normalise)
  
  
}