#' Function to plot output from ImmuneLENS
#'
#' @param vdj.region.df data frame containing coverage values by position
#' @param vdj.gene VDJ gene to be plotted
#' @param hg19_or_38 hg19 or hg38 version of genome
#' @param GC_correct TRUE or FALSE, whether to correct for GC content
#' @param GC_mode GC correct running 'simultaneous' to linear model or 'prior'
#' @param sample_name name of sample run
#' @param median.k rolling median window
#' @param median.thresh threshold to remove exons with low coverage
#' @param element.t.size size of text
#' @param point.size size of text
#' @param num.points Number of data points
#' @param customSeg Custom df defining normalisation segments
#' @param customFasta Custom fasta file
#' @param customFlag Custom list of flagged exons to remove
#' @param restriktor.absval parameter for fitting model with restriktor
#' @param ylims limits for y axis
#' @param removed_flag Whether to removed flagged WGS locations (default = TRUE)
#' @param output_df If true return data frame used for plotting instead of plot
#' @param classSwitch Whether to include class switching in model for IGH (default = TRUE)
#' @param customNorm Custom dataframe with average values to normalise for, summarised over 100bp (default = NULL)
#' @return data frame of TCRA T cell fractions with 95% CI
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth geom_segment geom_vline annotate theme_bw theme xlab ylab ylim scale_colour_manual ggtitle element_text
#' @importFrom tidyr pivot_longer
#' @importFrom ggnewscale new_scale_color
#' @name plotImmuneLENS
#' @export

plotImmuneLENS <- function(vdj.region.df, vdj.gene = 'TCRA',
                                 hg19_or_38 = 'hg38',GC_correct = TRUE,
                                 GC_mode = 'simultaneous',
                                 median.k = 50, median.thresh = 15,
                                 sample_name = 'test',
                                 num.points = 3000,
                                 ylims = NULL,element.t.size = NULL,
                                 point.size = 0.5,
                                 customSeg = NULL,
                                 customFasta = NULL,
                                 customFlag = NULL,
                                 restriktor.absval = 1e-5,
                                 removed_flag = TRUE,
                                 output_df = FALSE,
                                 classSwitch = TRUE,
                                 customNorm = NULL){

  # Checks on input data
  if(!vdj.gene %in% c('TCRA','TCRB','TCRG','IGH')){
    stop('Plotting currently only supported for TCRA, TCRB, TCRG or IGH')
  }
  if(!hg19_or_38 %in% c('hg19','hg38')){
    stop('hg19_or_38 must be either "hg19" or "hg38"')
  }
  if(!is.logical(GC_correct)){
    stop('GC.correct must be logical')
  }
  if(dim(vdj.region.df)[2] != 2 | !all(colnames(vdj.region.df) %in% c('pos','reads'))){
    stop('Expected vdf.region.df to be a data frame of 2 columns called pos and reads, is this not the output of loadCov()?')
  }

  # Sort out visible binding issue
  GC.correct <- segName <- start <- end2 <- hgnc_symbol <- NULL
  start_position_hg38 <- end_position_hg38 <- start_position_hg19 <- end_position_hg19 <- NULL
  X2 <- segType <- end <- prev_segType <- next_segType <- start2 <- NULL
  segment <- logR <- position <- logRgroup <- next_position <- pos <-  NULL
  Ratio.gc.correct <- next_logR <- exon2 <-  NULL
  GC_model_value <- Ratio <- exon.gc <- exon.gc2 <- smooth.gc <- smooth.gc2 <- NULL


  cov_df <- vdj.region.df

  if(!is.null(customFasta)){
    VDJ_fasta <- seqinr::read.fasta(customFasta)
  }else{
    if(hg19_or_38 == 'hg19'){
      # data("all_fasta")
      VDJ_fasta <- all_fasta[[vdj.gene]]
    }else{
      # data("all_fasta_hg38")
      VDJ_fasta <- all_fasta_hg38[[vdj.gene]]
    }
  }



  exclude.segs <- c('TRDD1','TRDD2','TRDD3','TRDC','TRAV11','TRAV8_4','TRAV8_5',
                    'TRAV15','TRAV31','TRAV7','TRAV9_1','TRAV18','TRAV8_7',
                    'TRAV28','TRAV32','TRAV33','TRAV37','TRDV2','TRDJ1','TRDJ2',
                    'TRDJ3','TRDJ4','TRAJ55','TRAJ51','TRAJ2','TRAJ1','TRAJ60',
                    'TRAJ59','TRAJ24','TRAJ14','TRBV7_1','TRBV6_2','TRBV8_1',
                    'TRBV5_2','TRBV8_2','TRBV12_3','TRBVA','TRBV26','TRBVB',
                    'TRBV1','TRBV5_3','TRB10_1','TRBV7_4','TRBV7_5','TRBV6_7',
                    'TRBV7_7','TRBV6_8','TRBV5_7','TRBV17','TRBV22_1','TRBC2',
                    'TRBV30','TRBD1','TRBC1','TRGC2','TRGC1',
                    'IGHGP','IGHEP1','IGHD',
                    # 'IGHA2','IGHE',
                    # 'IGHG4','IGHG2','IGHA1','IGHG1','IGHG3',
                    # 'IGHD','IGHM',
                    'IGHD1_26', 'IGHD6_25', 'IGHD5_24','IGHD4_23',
                    'IGHD3_22','IGHD2_21', 'IGHD1_20', 'IGHD6_19',
                    'IGHD5_18', 'IGHD4_17','IGHD3_16', 'IGHD2_15', 'IGHD1_14',
                    'IGHD6_13', 'IGHD5_12','IGHD4_11', 'IGHD3_10', 'IGHD3_9',
                    'IGHD2_8', 'IGHD1_7','IGHD6_6' , 'IGHD5_5', 'IGHD4_4',
                    'IGHD3_3', 'IGHD2_2', 'IGHD1_1','IGHD7_27','IGLC7','IGLL5',
                    'IGLC1','IGLC2','IGLC3','IGLC4','IGLC5','IGLC6','IGKC')

  if(!classSwitch){
    exclude.segs <- c(exclude.segs,
                      c('IGHA2','IGHE','IGHG4','IGHG2','IGHA1','IGHG1','IGHG3',
                        'IGHD','IGHM'))
  }

  # 1. Run exonsImmuneLENS for GC corrected positions
  # Flagged positions have already been removed if set to TRUE in function
  vdj.example.df <- exonsImmuneLENS(cov_df,vdj.gene,
                                                   hg19_or_38,
                                                   GC_correct = GC_correct,
                                                   removed_flag = removed_flag,
                                                   customFlag = customFlag,
                                                   median.thresh = median.thresh,
                                                   median.k = median.k,
                                                   sample_name = sample_name,
                                                   customSeg = customSeg,
                                                   customFasta = customFasta)


  # 2. Get segment data frames based on model
  if(vdj.gene == 'TCRA'){
    if(hg19_or_38 == 'hg38'){
      # start = end_position, end = next.start
      segment.ranges <- tcra_seg_hg38_vdj[-c(1:4),] %>%
        dplyr::filter(!segName %in% exclude.segs) %>%
        dplyr::mutate(end2 = dplyr::lead(start) - 1) %>%
        dplyr::select(segName, start, end = end2) %>%
        dplyr::mutate(segName = as.character(segName))
      segment.ranges$end[dim(segment.ranges)[1]] <- max(tcra_seg_hg38_vdj$end[-seq(4)])
    }else{
      if(hg19_or_38 == 'hg19'){
        # start = end_position, end = next.start
        segment.ranges <- tcra_seg_hg19_vdj[-c(1:4),] %>%
          dplyr::filter(!segName %in% exclude.segs) %>%
          dplyr::mutate(end2 = dplyr::lead(start) - 1) %>%
          dplyr::select(segName, start, end = end2) %>%
          dplyr::mutate(segName = as.character(segName))
        segment.ranges$end[dim(segment.ranges)[1]] <- max(tcra_seg_hg19_vdj$end[-seq(4)])
      }else{
        stop('hg19_or_38 must be hg19 or hg38')
      }
    }
  }else{
    if(hg19_or_38 == 'hg38'){
      segment.ranges <- vdj.segments.list[[vdj.gene]] %>%
        dplyr::select(segName = hgnc_symbol, start = start_position_hg38,
                      end = end_position_hg38) %>%
        dplyr::filter(!is.na(start)) %>%
        dplyr::filter(!segName %in% exclude.segs) %>%
        dplyr::mutate(end2 = dplyr::lead(start) - 1) %>%
        dplyr::select(segName, start, end = end2) %>%
        dplyr::mutate(segName = as.character(segName))
      segment.ranges$end[dim(segment.ranges)[1]] <- max(vdj.segments.list[[vdj.gene]]$end_position_hg38, na.rm = TRUE)
    }else{
      if(hg19_or_38 == 'hg19'){
        segment.ranges <- vdj.segments.list[[vdj.gene]] %>%
          dplyr::select(segName = hgnc_symbol, start = start_position_hg19,
                        end = end_position_hg19) %>%
          dplyr::filter(!is.na(start)) %>%
          dplyr::filter(!segName %in% exclude.segs) %>%
          dplyr::mutate(end2 = dplyr::lead(start) - 1) %>%
          dplyr::select(segName, start, end = end2) %>%
          dplyr::mutate(segName = as.character(segName))
        segment.ranges$end[dim(segment.ranges)[1]] <- max(vdj.segments.list[[vdj.gene]]$end_position_hg19, na.rm = TRUE)
      }else{
        stop('hg19_or_38 must be hg19 or hg38')
      }
      }
  }

  # Add in code for other VDJ genes

  # Preparation code for getVDJfraction_WGS_segmodel
  vdj.chr.df <- data.frame(gene = c('TCRA','TCRB','TCRG','IGH','IGL','IGK', 'TCRD'),
                           chr = c('chr14','chr7','chr7','chr14','chr22','chr2','chr14'))

  if(is.null(customSeg)){
    # data("vdj_seg_list")
    seg.name <- paste0(vdj.gene, '_', hg19_or_38)
    vdj.seg <- vdj_seg_list[[seg.name]]
  }else{
    vdj.seg <- customSeg
  }

  if(!is.null(customFasta)){
    VDJ_fasta <- seqinr::read.fasta(customFasta)
  }else{
    if(hg19_or_38 == 'hg19'){
      # data("all_fasta")
      VDJ_fasta <- all_fasta[[vdj.gene]]
    }else{
      # data("all_fasta_hg38")
      VDJ_fasta <- all_fasta_hg38[[vdj.gene]]
    }
  }

  vdj.chr <- vdj.chr.df$chr[which(vdj.chr.df$gene == vdj.gene)]
  vdj.start <- vdj.seg[1,2]
  vdj.end <- vdj.seg[1,3]

  TCRA.exons <- data.frame(X1 = vdj.chr,
                           X2 = seq(from = vdj.start,
                                    to = vdj.end - 1000, by = 1000)) %>%
    dplyr::mutate(X3 = X2 + 999)
  vdj.end.but1 <- TCRA.exons$X3[dim(TCRA.exons)[1]] + 1

  TCRA.exons <- rbind(TCRA.exons,
                      data.frame(X1 = vdj.chr, X2 = vdj.end.but1, X3 = vdj.end))

  exons.selected <- TCRA.exons

  # Get GC content for exons
  if(hg19_or_38 == 'hg19' & vdj.gene == 'TCRA'){
    exon.adjust.loc <- 21999999
  }else{
    exon.adjust.loc <- vdj.start - 1
  }

  TCRA.exons.loc <- list()
  for(i in seq_len(dim(exons.selected)[1])){
    TCRA.exons.loc[[i]] <- c(exons.selected$X2[i] - exon.adjust.loc,
                             exons.selected$X3[i] - exon.adjust.loc)

  }
  # These are the GC content within the exons



  exons.gc.content <- exonwindowplot2(TCRA.exons.loc, VDJ_fasta[[1]],0)


  if(GC_correct == TRUE & GC_mode == 'prior'){

    # Make GC correct ratio the input
    vdj.example.df$origRatio <- vdj.example.df$Ratio
    vdj.example.df$Ratio <- vdj.example.df$Ratio.gc.correct

    vdj.example.fit.df <- getVDJfraction_ImmuneLENS(vdj.example.df[,c(1,2)],
                                                      vdj.gene, sample_name,
                                                      hg19_or_38,
                                                      GC.correct = FALSE,
                                                      restriktor.absval = restriktor.absval,
                                                      exons = exons.selected,
                                                      exonList = exons.gc.content,
                                                      gene.fasta = VDJ_fasta,
                                                      gene.fasta.start = exon.adjust.loc,
                                                      classSwitch = classSwitch,
                                                      customNorm = customNorm)

  }else{

    vdj.example.fit.df <- getVDJfraction_ImmuneLENS(vdj.example.df[,c(1,2)],
                                                      vdj.gene, sample_name,
                                                      hg19_or_38,
                                                      GC.correct = GC_correct,
                                                      restriktor.absval = restriktor.absval,
                                                      exons = exons.selected,
                                                      exonList = exons.gc.content,
                                                      gene.fasta = VDJ_fasta,
                                                      gene.fasta.start = exon.adjust.loc,
                                                      classSwitch = classSwitch,
                                                      customNorm = customNorm)

  }





  # 3. Prepare ranges object
  if(vdj.gene == 'TCRA'){
    if(hg19_or_38 == 'hg38'){
      tcra.ranges <- tcra_seg_hg38_vdj[-c(1:4),]
    }else{
      if(hg19_or_38 == 'hg19'){
        tcra.ranges <- tcra_seg_hg19_vdj[-c(1:4),]
      }
    }
    tcra.ranges <- tcra.ranges %>%
      dplyr::mutate(segName = gsub('-','_',segName)) %>%
      dplyr::filter(!segName %in% exclude.segs) %>%
      dplyr::mutate(segType = ifelse(grepl('V',segName),'V','J')) %>%
      dplyr::mutate(next_segType = ifelse(grepl('V',dplyr::lead(segName)),'V','J')) %>%
      dplyr::mutate(prev_segType = ifelse(grepl('V',dplyr::lag(segName)),'V','J')) %>%
      dplyr::mutate(start2 = ifelse(segType == 'V', end + 1,
                                    ifelse(prev_segType == 'V',
                                           dplyr::lag(end) + 1,
                                           dplyr::lag(start)))) %>%
      dplyr::mutate(end2 = ifelse(segType == 'V',
                                  ifelse(next_segType == 'J',
                                         dplyr::lead(start) - 1,
                                         dplyr::lead(end)),
                                  start - 1))
    next_j_seg <- which(tcra.ranges$next_segType == 'J')
    tcra.ranges$segName[next_j_seg[1]] <- paste0(tcra.ranges$segName[next_j_seg[1]], '_',
                                                 tcra.ranges$segName[next_j_seg[1]+1])

    tcra.ranges <- tcra.ranges[-c(next_j_seg[1]+1),]
    tcra.ranges <- tcra.ranges %>%
      dplyr::select(segName, start = start2, end = end2) %>%
      dplyr::mutate(segName = as.character(segName))

    }else{
    if(hg19_or_38 == 'hg38'){
      tcra.ranges <- vdj.segments.list[[vdj.gene]] %>%
        dplyr::select(segName = hgnc_symbol, start = start_position_hg38, end = end_position_hg38) %>%
        dplyr::filter(!is.na(start))
    }
    if(hg19_or_38 == 'hg19'){
      tcra.ranges <- vdj.segments.list[[vdj.gene]] %>%
        dplyr::select(segName = hgnc_symbol, start = start_position_hg19, end = end_position_hg19) %>%
        dplyr::filter(!is.na(start)) %>%
        dplyr::arrange(start)
    }

    tcra.ranges <- tcra.ranges %>%
      dplyr::mutate(segName = gsub('-','_',segName)) %>%
      dplyr::filter(!segName %in% exclude.segs) %>%
      dplyr::mutate(segType = ifelse(grepl('V',segName),'V','J')) %>%
      dplyr::mutate(next_segType = ifelse(grepl('V',dplyr::lead(segName)),'V','J')) %>%
      dplyr::mutate(prev_segType = ifelse(grepl('V',dplyr::lag(segName)),'V','J'))

    if(vdj.gene %in% c('TCRB')){
      tcra.ranges <- tcra.ranges %>%
        dplyr::mutate(start2 = ifelse(segType == 'V', end + 1,
                                      ifelse(prev_segType == 'V',
                                             dplyr::lag(end) + 1,
                                             dplyr::lag(start)))) %>%
        dplyr::mutate(end2 = ifelse(segType == 'V',
                                    ifelse(next_segType == 'J',
                                           dplyr::lead(start) - 1,
                                           dplyr::lead(end)),
                                    start - 1))
      next_j_seg <- which(tcra.ranges$next_segType == 'J')
      tcra.ranges$segName[next_j_seg[1]] <- paste0(tcra.ranges$segName[next_j_seg[1]], '_',
                                                   tcra.ranges$segName[next_j_seg[1]+1])

      tcra.ranges <- tcra.ranges[-c(next_j_seg[1]+1),]
      tcra.ranges <- tcra.ranges %>%
        dplyr::select(segName, start = start2, end = end2) %>%
        dplyr::mutate(segName = as.character(segName))
    }
    if(vdj.gene %in% c('TCRG','IGH')){
      tcra.ranges <- tcra.ranges %>%
        dplyr::mutate(start2 = ifelse(segType == 'J', end + 1,
                                      ifelse(prev_segType == 'J',
                                             dplyr::lag(end) + 1,
                                             dplyr::lag(start)))) %>%
        dplyr::mutate(end2 = ifelse(segType == 'J',
                                    ifelse(next_segType == 'V',
                                           dplyr::lead(start) - 1,
                                           dplyr::lead(end)),
                                    start - 1))
      next_v_seg <- which(tcra.ranges$next_segType == 'V')
      tcra.ranges$segName[next_v_seg[1]] <- paste0(tcra.ranges$segName[next_v_seg[1]], '_',
                                                   tcra.ranges$segName[next_v_seg[1]+1])

      tcra.ranges <- tcra.ranges[-c(next_v_seg[1]+1),]
      tcra.ranges <- tcra.ranges %>%
        dplyr::select(segName, start = start2, end = end2) %>%
        dplyr::mutate(segName = as.character(segName))
    }
  }



  # 4. Prepare to plot

  exclude.segs2 <- c('exon.gc','exon.gc2','smooth.gc','smooth.gc2')
  if(!is.null(customNorm)){
    exclude.segs2 <- c(exclude.segs2, setdiff(colnames(customNorm),'exon2'))
  }

  tcra.seg.sol.df <- vdj.example.fit.df %>%
    dplyr::filter(!segment %in% exclude.segs2) %>%
    dplyr::select(segName = segment, logR) %>%
    dplyr::left_join(tcra.ranges, 'segName') %>%
    dplyr::arrange(start)

  # For class switching:
  if(vdj.gene == 'IGH' & classSwitch){
    IGH.M <- data.frame(segName = 'IGHM',logR = 0,
                        start = 105856219, end = 105863258)
    tcra.seg.sol.df <- rbind(tcra.seg.sol.df,
                             IGH.M) %>%
      dplyr::arrange(start)
  }

  tcra.seg.sol.df <- rbind(rbind(data.frame(segName = 'start', logR = 0,
                                            start = vdj.seg$start[1], end = (tcra.ranges$start)[1] - 1),
                                 tcra.seg.sol.df),
                           data.frame(segName = 'end',logR = 0,
                                      start = tcra.ranges$end[length(tcra.ranges$end)] + 1,
                                      end =vdj.seg$end[1]))


  tcra.seg.sol.df <- tcra.seg.sol.df[,c('logR','start','end')] %>%
    tidyr::pivot_longer(cols = c('start','end'), values_to = c('position')) %>%
    dplyr::select(position, logR)


  tcra.seg.sol.df <- tcra.seg.sol.df %>%
    dplyr::mutate(logRgroup= ifelse(abs(logR - dplyr::lag(logR,default = 0)) > 0.001,1,0)) %>%
    dplyr::mutate(logRgroup = cumsum(logRgroup) %% 2) %>%
    dplyr::mutate(logRgroup = ifelse(logRgroup == 0,'group1','group2')) %>%
    dplyr::mutate(next_position = dplyr::lead(position),
                  next_logR = dplyr::lead(logR))

  getColGroup <- function(x){
    groupID = tcra.seg.sol.df %>%
      dplyr::filter(position <= x & next_position > x)  %>%
      dplyr::select(logRgroup) %>% `[[`(1)
    # groupID = ifelse(groupID == 'group1','group2','group1')
    return(groupID)
  }
  getColGroup_v <- Vectorize(getColGroup)

  output_summary <- summariseWGS(vdj.example.fit.df,vdj.gene = vdj.gene)[[1]]
  output.frac <- output_summary[1,grep('fraction',colnames(output_summary))]
  output.div <- output_summary[1,grep('shannon',colnames(output_summary))[1]]

  t_or_b <- ifelse(vdj.gene %in% c('TCRA','TCRB','TCRG'), 'T', 'B')
  vseg_name <- ifelse(vdj.gene %in% c('TCRA','TCRB','TCRG'),
                      paste0(gsub('TC','T',vdj.gene),'V'), paste0(vdj.gene,'V'))
  # Add in parameters from function for adjusting sizes
  # Add in text annontation (cell fraction + diversity + sample name?)

  # If GC_mode = 'prior' can use Ratio.gc.correct
  # If GC_mode = 'simultaneous' we need to use the solutions to the model
  # to calculate a correction value
  if(GC_correct == TRUE & GC_mode == 'simultaneous'){

    if(!is.null(customNorm)){
      norm.names <- setdiff(colnames(customNorm),'exon2')
      vdj.gc.sol <- vdj.example.fit.df %>%
        dplyr::filter(segment %in% c('exon.gc','exon.gc2','smooth.gc','smooth.gc2',
                                     norm.names))
    }else{
      vdj.gc.sol <- vdj.example.fit.df %>%
        dplyr::filter(segment %in% c('exon.gc','exon.gc2','smooth.gc','smooth.gc2'))
    }



    solution.col <- which(colnames(vdj.example.fit.df) == 'logR')
    exon.gc.value <- vdj.gc.sol[vdj.gc.sol$segment == 'exon.gc',solution.col]
    exon.gc2.value <- vdj.gc.sol[vdj.gc.sol$segment == 'exon.gc2',solution.col]
    smooth.gc.value <- vdj.gc.sol[vdj.gc.sol$segment == 'smooth.gc',solution.col]
    smooth.gc2.value <- vdj.gc.sol[vdj.gc.sol$segment == 'smooth.gc2',solution.col]


    # Create new adjusted value
    vdj.example.df <- vdj.example.df %>%
      dplyr::mutate(GC_model_value = exon.gc.value*exon.gc + exon.gc2.value*exon.gc2 +
               smooth.gc.value*smooth.gc + smooth.gc2.value*smooth.gc2)

    if(!is.null(customNorm)){
      norm.names <- setdiff(colnames(customNorm),'exon2')
      norm.values <- sapply(norm.names, function(x){
        vdj.gc.sol[vdj.gc.sol$segment == x,solution.col]
      })

      vdj.example.df <- vdj.example.df %>%
        mutate(exon2 = floor((pos - (vdj.start - 1))/100) +1 ) %>%
        left_join(customNorm, 'exon2')

      for(i in seq_len(length(norm.names))){
        sym_name <- rlang::sym(norm.names[i])
        vdj.example.df <- vdj.example.df %>%
          dplyr::mutate(GC_model_value = GC_model_value + norm.values[i]*!!sym_name)
      }

    }
    vdj.example.df <- vdj.example.df %>%
      dplyr::mutate(Ratio.gc.correct = Ratio - GC_model_value)
  }

  if(GC_correct == TRUE){
    p1 <- vdj.example.df[sample(seq(dim(vdj.example.df)[1]),num.points),] %>%
      dplyr::mutate(logRgroup = getColGroup_v(pos)) %>%
      ggplot(aes(pos, Ratio.gc.correct)) +
      geom_point(aes(col = logRgroup),alpha = 0.1, size = point.size) +
      scale_colour_manual(values = c('red', 'blue')) +
      xlab('Genome position') + ylab('Log Ratio') +
      ggnewscale::new_scale_color()+
      geom_segment(aes(x = position, y = logR, xend = next_position,
                       yend = next_logR, col = logRgroup),data = tcra.seg.sol.df,
                   size = 1) +
      scale_colour_manual(values = c('darkred', 'darkblue')) +
      theme_bw() +
      theme(legend.position = 'none',
      ) +
      ggtitle(paste0(sample_name,': ',vdj.gene,' ',t_or_b,' cell fraction = ',
                     sprintf('%.3f',output.frac),'\n','Shannon ', vseg_name,
                     ' diversity = ',
                     sprintf('%.3f',output.div)))
  }else{
    p1 <- vdj.example.df[sample(seq(dim(vdj.example.df)[1]),num.points),] %>%
      dplyr::mutate(logRgroup = getColGroup_v(pos)) %>%
      ggplot(aes(pos, Ratio)) +
      geom_point(aes(col = logRgroup),alpha = 0.1, size = point.size) +
      scale_colour_manual(values = c('red', 'blue')) +
      xlab('Genome position') + ylab('Log Ratio') +
      ggnewscale::new_scale_color()+
      geom_segment(aes(x = position, y = logR, xend = next_position,
                       yend = next_logR, col = logRgroup),data = tcra.seg.sol.df,
                   size = 1) +
      scale_colour_manual(values = c('darkred', 'darkblue')) +
      theme_bw() +
      theme(legend.position = 'none',
      ) +
      ggtitle(paste0(sample_name,': ',vdj.gene,' ',t_or_b,' cell fraction = ',
                     sprintf('%.3f',output.frac),'\n','Shannon ', vseg_name,
                     ' diversity = ',
                     sprintf('%.3f',output.div)))
  }


  if(!is.null(ylims)){
    p1 <- p1 + ylim(ylims)
  }
  if(!is.null(element.t.size)){
    p1 <- p1 + theme(text = element_text(size = element.t.size))
  }
    if(output_df){
      return(list(vdj.example.df,tcra.seg.sol.df))
    }else{
      print(p1)
      return(p1)
    }

}
