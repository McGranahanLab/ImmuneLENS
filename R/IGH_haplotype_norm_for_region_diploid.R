#' Calculate predicted IGH germline haplotype assuming possible high B cell content
#'
#' @param input_cov data frame of positions and coverage values for IGH loci
#' @param  region_df threshold for length of genomic regions used in normalisation method (default = 5kb)
#' @return  IGH coverage data frame with normalised values by predicted IGH haplotype
#' @name IGH_haplotype_norm_for_region_diploid 

IGH_haplotype_norm_for_region_diploid <- function(input_cov, region_df){

  start_kb <- end_kb <- reads_seg <- pos <- value <- reads <- update_reads <- NULL
  
  if(dim(region_df)[1] == 0){
    warning('No regions to normalise')
    return(input_cov)
  }
  
  bin_size <- region_df$bin_size[1]
  # Region has already been defined
  ranges_to_normalise_long <- region_df %>%
    dplyr::rowwise() %>%
    dplyr::transmute(
      pos = list(seq(start_kb * bin_size , end_kb * bin_size + (bin_size - 1))),
      value = reads_seg
    ) %>%
    tidyr::unnest(cols = c(pos))
  
  # If kb_bin in ranges to normalise update based on reads_seg_average_roll value
  
  update_cov <- input_cov %>% 
    dplyr::mutate(kb_bin = pos %/% bin_size ) %>%
    dplyr::left_join(ranges_to_normalise_long, 'pos') %>%
    dplyr::mutate(value = ifelse(is.na(value),1, value)) %>%
    dplyr::filter(value != 0) %>%
    dplyr::mutate(update_reads = reads*(1/value)) %>%
    dplyr::select(pos, reads = update_reads)
  
  return(update_cov)
  
}
