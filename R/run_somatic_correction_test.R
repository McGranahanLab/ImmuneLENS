#' QC test for predicted somatic copy number alterations in IGH loci (max logR change > 0.1 and kb change length > 100)
#'
#' @param somatic_regions_df data frame of predicted somatic copy number alterations in IGH loci
#' @return TRUE or FALSE based on QC test
#' @name run_somatic_correction_test

run_somatic_correction_test <- function(somatic_regions_df){
  if(dim(somatic_regions_df)[[1]] == 1){
    return(FALSE)
  }
  focal_start <- vdj_seg_list[['IGH_hg38']][2,2] %/% 1000
  focal_end <- vdj_seg_list[['IGH_hg38']][2,3] %/% 1000
  regions_without_focal <- somatic_regions_df %>% 
    # filter(!(start_kb <= focal_start & end_kb >= focal_end)) %>%
    filter(germline_reads > 0.1)
  
  kb_change_length <- regions_without_focal  %>% filter(reads_seg != 1) %>% 
    select(kb_length) %>% `[[`(1) %>% sum()
  
  max_logR_change <- max(abs(1-regions_without_focal$logR),na.rm = TRUE)
  #range_logR <- max_logR - min(regions_without_focal$logR, na.rm = TRUE)
  
  if(max_logR_change> 0.1 & kb_change_length > 100){
    return(TRUE)
  }else{
    return(FALSE)
  }
}
