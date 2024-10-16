#' Calculate median values of coverage within exons in rolling windows
#'
#' @param vdj.region.df data frame of positions and coverage values
#' @param exons.selected Locations of exons
#' @param median.k rolling median window
#' @param median.thresh threshold to remove exons with low coverage
#' @param exons.to.use option to manually select exons
#' @return Data frame of positions and rolling median of coverage values within chosen exons
#' @importFrom stats median
#' @name medianExonCoverage

medianExonCoverage <- function(vdj.region.df, exons.selected, median.k = 50, median.thresh = 15, exons.to.use){

  pos <- NULL
  exon_id <- start <- end <- reads <-  NULL
  

  # Filter for positions within exons as expected and apply median filter and remove any exons required
  # Convert to data.table for faster processing
  vdj.region.dt <- data.table::as.data.table(vdj.region.df)
  exons.selected.dt <- data.table::as.data.table(exons.selected)
  
  # Prepare data for optimised approach
  vdj.region.dt$start <- vdj.region.dt$pos
  vdj.region.dt$end <- vdj.region.dt$pos
  # vdj.region.dt[, `:=`(start = pos, end = pos)]
  
  exons_subset <- exons.to.use
  exons_subset.dt <- exons.selected.dt[exons_subset,]
  data.table::setnames(exons_subset.dt, old = c("X2", "X3"), new = c("start", "end"))
  exons_subset.dt$exon_id <- seq_len(dim(exons_subset.dt)[1])
  data.table::setkey(vdj.region.dt, start, end)
  data.table::setkey(exons_subset.dt, start, end)
  
  overlaps <- data.table::foverlaps(vdj.region.dt, exons_subset.dt, type = "within", nomatch = 0L)
  reads_dt <- overlaps[, .(pos = list(pos), reads = list(reads)), by = exon_id]
  
  vdj.region.df.filt.exons.median <- lapply(seq_len(dim(reads_dt)[1]), function(x){
    tmp <- data.frame(pos = reads_dt$pos[[x]], reads = reads_dt$reads[[x]])
    tmp$reads <- medianFilter(tmp$reads, median.k)
    return(tmp)})
  
  # Remove any exons with very low values (suspected exon failure - possible 100% T cell)
  # Threshold of < 15 may not be appropriate on low coverage data sets/genomic regions
  median.values.exons <- sapply(vdj.region.df.filt.exons.median,
                                function(x) median(x$reads, na.rm = TRUE))
  exon.remove <- which(median.values.exons < median.thresh | is.na(median.values.exons))
  if(length(exon.remove) > 0){
    vdj.region.df.filt.exons.median <- vdj.region.df.filt.exons.median[-exon.remove]
  }
  vdj.region.df.filt.exons.median <- data.frame(data.table::rbindlist(vdj.region.df.filt.exons.median))
  
  # Get rid of any repeated rows that might exist - in theory I don't think they will
  vdj.region.df.filt.exons.median <- dplyr::distinct(vdj.region.df.filt.exons.median)
  return(list(vdj.region.df.filt.exons.median, exon.remove))
  
}
