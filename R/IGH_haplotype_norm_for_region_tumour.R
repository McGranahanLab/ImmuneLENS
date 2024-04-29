#' Function to normalise IGH coverage values by predicted IGH germline haplotype regions for a tumour sample
#'
#' @param input_cov data frame of positions and coverage values for IGH loci
#' @param  region_df predicted IGH germline haplotype regions
#' @param norm_type type of normalisation to use (default = 'benchmark')
#' @param gc_correct logical, should the normalisation be GC corrected? (default = FALSE)
#' @return  IGH coverage data frame with normalised values by predicted IGH germline haplotype applied to a tumour sample
#' @name IGH_haplotype_norm_for_region_diploid 

IGH_haplotype_norm_for_region_tumour <- function(input_cov, region_df,
                                                 norm_type = 'benchmark', 
                                                 gc_correct = FALSE){
  reads_seg <- start_kb <- end_kb <- pos <- value <- reads <- update_reads <- NULL
  
  # We do not assume any knowledge on the tumour purity and allele specific copy-number but for each 
  # region will calculate the median change in coverage and use that
  # as the basis of the normalisation values. 
  if(dim(region_df)[1] == 0){
    warning('No regions to normalise')
    return(input_cov)
  }
  
  bin_size <- region_df$bin_size[1]
  if(norm_type == 'benchmark'){
    if(gc_correct){
      input_cov_gc <- IGH_gc_correct(input_cov)
      norm_values <- sapply(seq_len(dim(region_df)[1]),
                            function(x) {
                              calc_norm_for_region(input_cov_gc, region_df$start_kb[x],
                                                   region_df$end_kb[x], bin_size = bin_size)})
    }else{
      norm_values <- sapply(seq_len(dim(region_df)[1]),
                            function(x) {
                              calc_norm_for_region(input_cov, region_df$start_kb[x],
                                                   region_df$end_kb[x], bin_size = bin_size)})
    }
    
  }
  if(norm_type == 'local'){
    if(gc_correct){
      input_cov_gc <- IGH_gc_correct(input_cov)
      norm_values <- sapply(seq_len(dim(region_df)[1]),
                            function(x) {
                              calc_norm_for_region_alt(input_cov_gc, region_df$start_kb[x],
                                                       region_df$end_kb[x], bin_size = bin_size,
                                                       all_regions = region_df)})
    }
    else{
      norm_values <- sapply(seq_len(dim(region_df)[1]),
                            function(x) {
                              calc_norm_for_region_alt(input_cov, region_df$start_kb[x],
                                                       region_df$end_kb[x], bin_size = bin_size,
                                                       all_regions = region_df)})
    }
    
  }
  
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
      pos = list(seq(start_kb * bin_size, end_kb * bin_size + (bin_size-1))),
      value = reads_seg,
      norm = norm
    ) %>%
    tidyr::unnest(cols = c(pos))
  
  
  # If kb_bin in ranges to normalise update based on reads_seg_average_roll value
  
  update_cov <- input_cov %>% 
    dplyr::mutate(kb_bin = pos %/% bin_size) %>%
    dplyr::left_join(ranges_to_normalise_long, 'pos') %>%
    dplyr::mutate(value = ifelse(is.na(value), 1,  value)) %>%
    dplyr::filter(value != 0) %>%
    dplyr::mutate(update_reads = ifelse(value == 1, reads, reads*norm)) %>%
    dplyr::select(pos, reads = update_reads)
  
  return(update_cov)
  
}

calc_norm_for_region <- function(input_cov, start_pos, end_pos,bin_size, kb_range = 2500){
  pos <- reads <- NULL
  # Median coverage in region:
  baseline_median <- input_cov %>%
    dplyr::filter(pos < 105588395 | pos > 106779812)  %>%
    dplyr::select(reads) %>% `[[`(1) %>% median(na.rm = TRUE)
  baseline_end <- input_cov %>%
    dplyr::filter( pos > 106779812)  %>%
    dplyr::select(reads) %>% `[[`(1) %>% median(na.rm = TRUE)
  baseline_start <- input_cov %>%
    dplyr::filter(pos < 105588395 )  %>%
    dplyr::select(reads) %>% `[[`(1) %>% median(na.rm = TRUE)
  baseline_all <- median(input_cov$reads, na.rm = TRUE)
  
  # Check - is start/end baseline different? This is a big problem!
  if(abs(baseline_start - baseline_end) < round(median(baseline_all, na.rm = TRUE)/10)){
    baseline_to_use <- baseline_median
  }else{
    # Which is the closest to the overall median?
    start_diff <- abs(baseline_all - baseline_start)
    end_diff <- abs(baseline_all - baseline_end)
    if(start_diff < end_diff){
      baseline_to_use <- baseline_start
    }else{
      baseline_to_use <- baseline_end
    }
  }
  
  start_pos <- start_pos*bin_size
  end_pos <- end_pos*bin_size + (bin_size - 1)
  region_median <- input_cov %>%
    dplyr::filter(pos >= start_pos) %>%
    dplyr::filter(pos <= end_pos) %>%
    dplyr::select(reads) %>% `[[`(1) %>% median(na.rm = TRUE)

  
  # Before vs region:
  region_change <- baseline_to_use/region_median

  norm.value = region_change
  
  return(norm.value)
  
}

calc_norm_for_region_alt <- function(input_cov, start_pos, end_pos,bin_size,all_regions,kb_range = 250000){
  pos <- reads <- NULL
  input_cov_baseline <- input_cov
  for(i in seq_len(dim(all_regions)[1])){
    tmp_start <- all_regions$start_kb[i]*bin_size
    tmp_end <- (all_regions$end_kb[i]*bin_size) + bin_size -1 
    input_cov_baseline <- input_cov_baseline  %>%
      dplyr::filter(pos < tmp_start | pos > tmp_end)
  }
  
  start_pos <- start_pos*bin_size
  end_pos <- end_pos*bin_size + (bin_size - 1)
  region_median <- input_cov %>%
    dplyr::filter(pos >= start_pos) %>%
    dplyr::filter(pos <= end_pos) %>%
    dplyr::select(reads) %>% `[[`(1) %>% median(na.rm = TRUE)
  # Median coverage in kb_range before region:
  pre_region_median <- input_cov_baseline %>%
    dplyr::filter(pos >= start_pos - kb_range) %>%
    dplyr::filter(pos <= start_pos) %>%
    dplyr::select(reads) %>% `[[`(1) 
  # 
  # # Median coverage in kb_range after region:
  post_region_median <- input_cov_baseline %>%
    dplyr::filter(pos >= end_pos + bin_size -1) %>%
    dplyr::filter(pos <= end_pos + bin_size + kb_range) %>%
    dplyr::select(reads) %>% `[[`(1)
  
  pre_post_region_median <- median(c(pre_region_median, post_region_median),na.rm = TRUE)
  
  region_change = (pre_post_region_median)/region_median
  
  norm.value = region_change
  
  return(norm.value)
  
}
