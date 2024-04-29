#' Calculate normalised coverage file for predicted IGH germline haplotype
#'
#' @param cov.df data frame of positions and coverage values for IGH loci
#' @param kb_len_threshold threshold for length of genomic regions used in normalisation method (default = 5kb)
#' @param bin_size bins used for normalisation (default = 1000)
#' @param highBcell threshold to remove exons with low coverage
#' @return List of 1. normalised coverage dataframe and 2. regions used for normalisation
#' @name IGH_haplotype_norm_fun

IGH_haplotype_norm_fun <- function(cov.df, kb_len_threshold = 5, bin_size = 1000, highBcell = FALSE) {
  cov.df.gc <- IGH_gc_correct_alt(cov.df)
  
  if(highBcell){
    regions_gl_df <- calc_haplotype_regions_de_novo_highBcell(cov.df.gc,
                                                              kb_len_threshold = kb_len_threshold,
                                                              bin_size = 1000)
  }else{
    regions_gl_df1 <- calc_haplotype_regions_de_novo_diploid(cov.df.gc,
                                                             kb_len_threshold = kb_len_threshold,
                                                             bin_size = 1000)
    
  }
  
  cov.df_update1 <- IGH_haplotype_norm_for_region_diploid(cov.df, regions_gl_df1)
  
  return(list(cov.df_update1, regions_gl_df1))
}