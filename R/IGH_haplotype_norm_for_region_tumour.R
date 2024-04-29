IGH_haplotype_norm_for_region_tumour <- function(input_cov, region_df,
                                                 norm_type = 'benchmark', 
                                                 gc_correct = FALSE){
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
      input_cov_gc <- IGH_gc_correct_alt(input_cov)
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
      input_cov_gc <- IGH_gc_correct_alt(input_cov)
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