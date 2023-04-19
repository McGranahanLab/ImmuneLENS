#' Adjust TCRA T cell fraction based on tumour purity and copy number
#'
#' @param out.df Output containing unadjusted T/B cell fraction
#' @param purity Tumour purity, between 0 and 1
#' @param local.cn Tumour copynumber at VDJ gene locus
#' @param vdj.gene VDJ gene of interest
#' @importFrom rlang :=
#' @return data frame of TCRA scores with purity/copy number adjusted outputs
#' @name adjustImmuneLENS
#' @export

adjustImmuneLENS <- function(out.df, purity, local.cn, vdj.gene){
  if(vdj.gene %in% c('TCRA','TCRB','TCRG')){
    fraction_col <- rlang::sym(paste0(vdj.gene,'.tcell.fraction'))
    fraction_col_adj <- rlang::sym(paste0(vdj.gene,'.tcell.fraction.adj'))
  }else if(vdj.gene == 'IGH'){
    fraction_col <- rlang::sym(paste0(vdj.gene,'.bcell.fraction'))
    fraction_col_adj <- rlang::sym(paste0(vdj.gene,'.bcell.fraction.adj'))
  }
  if(any(colnames(out.df) == 'segment.fraction')){
    fraction_col <- rlang::sym('segment.fraction')
    fraction_col_adj <- rlang::sym('segment.fraction.adj')
  }

  cn_col <- rlang::sym(paste0(vdj.gene,'.cn'))

  # solve visible binding issue
  rawRatio <- maxPossible <- highTcellFlag <- NULL


  # check purity is correct
  if(any(purity > 1) | any(purity < 0)) stop('purity needs to be between 0 and 1')

  out.df$purity <- purity
  out.df[[paste0(vdj.gene,'.cn')]] <- local.cn

  out.df <- out.df %>%
    dplyr::mutate(rawRatio = 1-!!fraction_col) %>%
    dplyr::mutate(!!fraction_col_adj := 1 - ((1-purity+(purity*!!cn_col)/2)*rawRatio) - purity + ((purity*!!cn_col)/2)) %>%
    dplyr::mutate(maxPossible = 1- purity) %>%
    dplyr::mutate(highTcellFlag = !!fraction_col_adj > maxPossible)

  # Remove rawRatio column
  out.df <- out.df[,-which(colnames(out.df) == 'rawRatio')]

  return(out.df)

}
