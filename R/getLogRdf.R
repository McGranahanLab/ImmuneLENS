#' Method for calculating log ratio from single region
#'
#' @param region.df coverage values with a specific region
#' @param segs Location of segments used for normalisation and focal region
#' @param minCov Minimum GC coverage required
#' @return Calculated single log ratio based on segments
#' @importFrom stats median
#' @name getLogRdf

getLogRdf <-  function(region.df, segs, minCov = 0){

  pos <- reads <- NULL

  col_input <- 'reads'
  col.sym <- rlang::sym(col_input)
  # For random locations with no VDJ effect use beginning and end of TCRA
  tumour.random.covs1 <- region.df %>%
    dplyr::filter(pos <= segs[3,]$end) %>%
    dplyr::filter(reads >=  minCov) %>% dplyr::select(!!col.sym) %>% `[[`(1)
  tumour.random.covs2 <- region.df %>%
    dplyr::filter(pos >= segs[4,]$start) %>%
    dplyr::filter(reads >=  minCov) %>% dplyr::select(!!col.sym) %>% `[[`(1)
  tumour.random.covs <- c(tumour.random.covs1, tumour.random.covs2)
  
  if(length(tumour.random.covs) == 0) stop(paste0('No positions with coverage in the normalisation regions: ',
                                                  segs$start[3], '-',segs$end[3],' and ',
                                                  segs$start[4], '-',segs$end[4]))
  
  # Use median of these values for normalisation to get "logR"
  n1 <- median(tumour.random.covs)

  # If median = 0 (low coverage cases - instead use mean)
  if(n1 == 0){
    n1 <- mean(tumour.random.covs)
    if(n1 == 0) stop('All normalisation positions have 0 coverage. Impossible to calculate logR')
  }
  
  
  # Select coverage values across TCRA
  tumour.test.covs1 <- region.df %>%
    dplyr::filter(pos >= segs[1,]$start & pos <= segs[1,]$end) %>%
    dplyr::filter(reads >=  minCov) %>% dplyr::select(!!col.sym) %>% `[[`(1)

  # Select positions
  tumour.test.covs.pos <- region.df %>%
    dplyr::filter(pos >= segs[1,]$start & pos <= segs[1,]$end) %>%
    dplyr::filter(reads >=  minCov) %>% dplyr::select(pos) %>% `[[`(1)

  # Adjust for bin width
  # tumour.test.covs.pos <- tumour.test.covs.pos + bin.width/2

  # Get LogR df
  tumour.logR <- data.frame(pos = tumour.test.covs.pos,
                            Ratio = log2(tumour.test.covs1/n1))

  return(tumour.logR)

}
