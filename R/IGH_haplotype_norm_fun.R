#' Calculate normalised coverage file for the predicted IGH germline haplotype
#'
#' @param cov.df data frame of positions and coverage values for IGH loci
#' @param kb_len_threshold threshold for length of genomic regions used in normalisation method (default = 5kb)
#' @param bin_size bins used for normalisation (default = 1000)
#' @param highBcell Is the sample high in B cells (> 0.25 fraction)? (default = FALSE)
#' @return List of 1. normalised coverage dataframe and 2. regions used for normalisation
#' @name IGH_haplotype_norm_fun
#' @export


IGH_haplotype_norm_fun <- function(cov.df, kb_len_threshold = 5, bin_size = 1000, highBcell = FALSE) {
  cov.df.gc <- IGH_gc_correct(cov.df)
  
  if(highBcell){
    regions_gl_df <- calc_haplotype_regions_de_novo_highBcell(cov.df.gc,
                                                              kb_len_threshold = kb_len_threshold)
  }else{
    regions_gl_df1 <- calc_haplotype_regions_de_novo_diploid(cov.df.gc,
                                                             kb_len_threshold = kb_len_threshold,
                                                             bin_size = 1000)
    
  }
  
  cov.df_update1 <- IGH_haplotype_norm_for_region_diploid(cov.df, regions_gl_df1)
  
  return(list(cov.df_update1, regions_gl_df1))
}