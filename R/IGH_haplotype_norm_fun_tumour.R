#' Calculate normalised coverage file for the predicted IGH germline haplotype
#'
#' @param cov.df data frame of positions and coverage values for IGH loci of matched germline sample
#' @param cov.df.tumour data frame of positions and coverage values for IGH loci of tumour sample
#' @param regions_gl_df1 data frame of regions used for normalisation from germline sample
#' @param gc_correct logical, should the normalisation be GC corrected? (default = TRUE)
#' @param run_breakpoint_check logical, should the normalisation be checked for breakpoints? (default = TRUE)
#' @param threshold_to_norm threshold for normalisation (default = 0)
#' @param norm_type_to_use type of normalisation to use (default = 'benchmark')
#' @return List of four items. 1. normalised coverage without somatic normalisation. 2. normalisation coverage with somatic normalisaton. 3. Somatic regions predicted and used for normalisation. 4. TRUE/FALSE of whether somatic regions passed QC
#' @name IGH_haplotype_norm_fun_tumour
#' @export

IGH_haplotype_norm_fun_tumour <- function(cov.df, cov.df.tumour, regions_gl_df1, ,
                                          gc_correct = TRUE,
                                          run_breakpoint_check = TRUE, threshold_to_norm = 0,
                                          norm_type_to_use = 'benchmark') {
  
  rpart_cp = 0.005
  
  somatic_regions_df <- somatic_cna_IGH_caller(cov.df, cov.df.tumour,
                                               gc_correct = gc_correct, 
                                               rpart_cp = rpart_cp, 
                                               run_breakpoint_check = run_breakpoint_check)
  
  somatic_regions_df$correct <- run_somatic_correction_test(somatic_regions_df)
  
  if (somatic_regions_df$correct[1]) {
    tumour_cov_somatic_update <- IGH_somatic_norm_for_region(cov.df.tumour, 
                                                             region_df = somatic_regions_df, 
                                                             threshold_to_norm = threshold_to_norm)
    cov.df.tumour2 <- tumour_cov_somatic_update
  } else {
    cov.df.tumour2 <- cov.df.tumour
  }
  
  cov.df_update2a <- IGH_haplotype_norm_for_region_tumour(cov.df.tumour, regions_gl_df1, 
                                                          norm_type = norm_type_to_use,
                                                          gc_correct = gc_correct)
  
  cov.df_update2b <- IGH_haplotype_norm_for_region_tumour(cov.df.tumour2, regions_gl_df1, 
                                                          norm_type = norm_type_to_use,
                                                          gc_correct = gc_correct)
  
  # Return: output without any somatic correction, output with somatic correction, somatic regions, somatic correction applied
  return(list(cov.df_update2a = cov.df_update2a, 
              cov.df_update2b = cov.df_update2b, 
              somatic_regions_df = somatic_regions_df, 
              somatic_correction_applied = somatic_regions_df$correct[1]))
}
