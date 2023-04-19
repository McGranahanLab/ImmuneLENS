#' Function to calculate distance in V or J segment usage between two samples
#'
#' @param out_df1 data frame output from runTcellExtrect_WGS containing segment fractions for sample 1
#' @param out_df2 data frame output from runTcellExtrect_WGS containing segment fractions for sample 2
#' @param V_or_J plot V or J segment proportions
#' @return JSD_dist
#' @name calcJSD_distance
#' @export
calcJSD_distance <- function(out_df1, out_df2, V_or_J){
  if(!V_or_J %in% c('V','J')){
    stop('V_or_J must be either "V" or "J"')
  }

  segment.prop <- segment.fraction <- prop1 <- prop2 <- n <- NULL
  segment <- JSD_1 <- JSD_2 <- NULL

  out_df1 <- out_df1 %>%
    dplyr::filter(grepl(V_or_J, segment)) %>%
    dplyr::mutate(segment.prop = (segment.fraction)/sum(segment.fraction))  %>%
    dplyr::mutate(segment.prop = segment.prop + 1e-5)%>%
    dplyr::select(segment, prop1 = segment.prop)

  out_df2 <- out_df2 %>%
      dplyr::filter(grepl(V_or_J, segment)) %>%
      dplyr::mutate(segment.prop = (segment.fraction)/sum(segment.fraction)) %>%
      dplyr::mutate(segment.prop = segment.prop + 1e-5)  %>%
      dplyr::select(segment, prop2 = segment.prop)

  out_df_all <- dplyr::inner_join(out_df1, out_df2, 'segment') %>%
    dplyr::mutate(n = 0.5 * (prop1 + prop2)) %>%
    dplyr::mutate(JSD_1 = prop1 * log(prop1/n)) %>%
    dplyr::mutate(JSD_2 = prop2 * log(prop2/n))

  JSD_dist <- 0.5*(sum(out_df_all$JSD_1) + sum(out_df_all$JSD_2))
  return(JSD_dist)

}
