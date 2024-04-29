#' Function to normalise IGH coverage values by predicted somatic copy number alterations
#'
#' @param input_cov data frame of coverage values in IGH loci
#' @param region_df data frame of predicted copy number of IGH loci regions for normalisations
#' @param threshold_to_norm threshold for genomic regions to be normalised
#' @return Updated coverage data frame normalised by predicted somatic copy number alterations
#' @name IGH_somatic_norm_for_region


IGH_somatic_norm_for_region <- function(input_cov, region_df, threshold_to_norm = 0.2){
  start_kb <- end_kb <- reads_seg <- pos <- value <- reads <- update_reads <- NULL
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
    dplyr::mutate(value = ifelse(abs(1-value) < threshold_to_norm, 1, value)) %>%
    dplyr::mutate(update_reads = reads*(1/value)) %>%
    dplyr::select(pos, reads = update_reads)
  
  return(update_cov)
  
}
