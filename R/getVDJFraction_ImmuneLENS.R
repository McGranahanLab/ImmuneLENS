#' Function to run ImmuneLENS for WGS for individual V and J segments
#'
#' @param test.logR data frame containing logR values to be fitted to constrained model
#' @param vdj.gene V(D)J gene to use
#' @param sample_name name of the sample
#' @param  hg19_or_38 version of the genome
#' @param GC.correct whether GC.correct is being used
#' @param restriktor.absval Parameter for use in Restriktor
#' @param exons Location of exons
#' @param exonList GC content of exons
#' @param gene.fasta FASTA file for VDJ gene e.g TCRA
#' @param gene.fasta.start offset for FASTA file, e.g 21999999 for TCRA in hg19
#' @param sliding number of bp for gc windows (default = 1000)
#' @param classSwitch Whether to include class switching in model for IGH (default = TRUE)
#' @param customNorm Custom dataframe with average values to normalise for, summarised over 100bp (default = NULL)
#' @return data frame of TCRA T cell fractions for VDJ segments
#' @name getVDJfraction_ImmuneLENS
#'
getVDJfraction_ImmuneLENS <- function(test.logR, vdj.gene, sample_name,
                                        hg19_or_38 = 'hg38', GC.correct = TRUE,
                                        restriktor.absval = 1e-5,
                                        exons = NULL,
                                        exonList = NULL,
                                        gene.fasta = NULL,
                                        gene.fasta.start = NULL,
                                        sliding = 1000,
                                        classSwitch = TRUE,
                                        customNorm = NULL){


  # Solve latest binding issues
  segName <- segType <- end <- prev_segType <- NULL
  start <- next_segType <- start2 <- end2 <-NULL
  start_position_hg19 <- end_position_hg19 <- NULL
  hgnc_symbol <- start_position_hg38 <- NULL
  end_position_hg38 <- Ratio <- exon.gc2 <- smooth.gc2 <- comparison <- NULL
  segName2 <- constraint <- logR <- segment <- VDJ.tcell.fraction <- NULL
  X2 <- NULL
  pos <- GC <- exon.gc <- smooth.gc <- NULL

  if(GC.correct){

    TCRA.gc.df <- slidingwindowplot_alt(sliding, gene.fasta[[1]])
    TCRA.gc.df$pos <- TCRA.gc.df$loc + gene.fasta.start
    gam.model <- mgcv::gam(GC~s(pos, bs = 'cs'), data = TCRA.gc.df)

    # This is the function to get the smoothed GC content at any position
    get_gc_prediction <- function(x){
      mgcv::predict.gam(gam.model,newdata = data.frame(pos = x))}

    test.logR <- test.logR %>%
      dplyr::mutate(exon = exonPosFun_v(pos, exons)) %>%
      dplyr::left_join(exonList, 'exon') %>%
      dplyr::rename(exon.gc = GC) %>%
      dplyr::mutate(exon.gc2 = exon.gc^2) %>%
      dplyr::mutate(smooth.gc = get_gc_prediction(pos))  %>%
      dplyr::mutate(smooth.gc2 = smooth.gc^2)


    # Add custom norm columns
    if(!is.null(customNorm)){

      start_pos_1 <- vdj_seg_list[[paste0(vdj.gene,'_',hg19_or_38)]]$start[1] - 1

      test.logR <- test.logR %>%
        mutate(exon2 = floor((pos - start_pos_1)/100) +1 ) %>%
        left_join(customNorm,'exon2')

      norm_col_names <- setdiff(colnames(customNorm),'exon2')

    }

  }


  create_col_seg <- function(input_df, start_pos, end_pos,col.name){
    input_df <- input_df %>%
      dplyr::mutate(X1= ifelse(pos >= start_pos & pos <= end_pos, 1, 0))
    colnames(input_df)[which(colnames(input_df) == 'X1')] <- col.name
    return(input_df)
  }

  # Vector of segments to exclude from model
  # e.g. pseudogenes or rarely used in human population

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

  if(vdj.gene == 'TCRA'){
    if(hg19_or_38 == 'hg38'){
      # data("tcra_seg_hg38_vdj_version")
      segment.ranges <- tcra_seg_hg38_vdj[-c(1:4),]
    }else{
      if(hg19_or_38 == 'hg19'){
        # data("tcra_seg_hg19_vdj_version")
        segment.ranges <- tcra_seg_hg19_vdj[-c(1:4),]
        }else{
          stop('hg19_or_38 must be hg19 or hg38')
        }
    }
    segment.ranges <- segment.ranges %>%
      dplyr::mutate(segName = gsub('-','_', segName)) %>%
      dplyr::filter(!segName %in% exclude.segs) %>%
      dplyr::mutate(segType = ifelse(grepl('V',segName),'V','J')) %>%
      dplyr::mutate(next_segType = ifelse(grepl('V',dplyr::lead(segName)),'V','J')) %>%
      dplyr::mutate(prev_segType = ifelse(grepl('V',dplyr::lag(segName)),'V','J')) %>%
      dplyr::mutate(start2 = ifelse(segType == 'V',end + 1,
                                    ifelse(prev_segType == 'V',
                                           dplyr::lag(end) + 1,dplyr::lag(start)))) %>%
      dplyr::mutate(end2 = ifelse(segType == 'V',
                                  ifelse(next_segType == 'J',
                                         dplyr::lead(start) - 1,dplyr::lead(end)),
                                  start -1))

    j_seg_loc <- which(segment.ranges$next_segType == 'J')
    new_seg_name <- paste0(segment.ranges$segName[j_seg_loc[1]],
                           '_', segment.ranges$segName[j_seg_loc[1] + 1])
    segment.ranges$segName[j_seg_loc[1]] <- new_seg_name
    segment.ranges <- segment.ranges[-c(j_seg_loc[1] + 1), ]
    segment.ranges <- segment.ranges %>%
      dplyr::select(segName, start = start2, end = end2) %>%
      dplyr::mutate(segName = as.character(segName))

  }else{
    if(hg19_or_38 == 'hg38'){
      if(vdj.gene %in% c('TCRG','IGH')){
        segment.ranges <- vdj.segments.list[[vdj.gene]] %>%
          dplyr::select(segName = hgnc_symbol, start = start_position_hg38,
                        end = end_position_hg38) %>%
          dplyr::mutate(segName = gsub('-','_', segName)) %>%
          dplyr::filter(!is.na(start)) %>%
          dplyr::filter(!segName %in% exclude.segs) %>%
          dplyr::mutate(segType = ifelse(grepl('V',segName),'V','J')) %>%
          dplyr::mutate(next_segType = ifelse(grepl('V',dplyr::lead(segName)),'V','J')) %>%
          dplyr::mutate(prev_segType = ifelse(grepl('V',dplyr::lag(segName)),'V','J')) %>%
          dplyr::mutate(start2 = ifelse(segType == 'J',end + 1,
                                        ifelse(prev_segType == 'J',
                                               dplyr::lag(end) + 1,dplyr::lag(start)))) %>%
          dplyr::mutate(end2 = ifelse(segType == 'J',
                                      ifelse(next_segType == 'V',
                                             dplyr::lead(start) - 1,dplyr::lead(end)),
                                      start -1))

        if(vdj.gene == 'IGH' & classSwitch){
          # Define class switching gene segments and update segment ranges
          class.switch.genes <- c('IGHA2','IGHE','IGHG4','IGHG2','IGHA1',
                                  'IGHG1','IGHG3','IGHM')
          segment.ranges <- segment.ranges %>%
            dplyr::mutate(segType = ifelse(segName %in% class.switch.genes,
                                    'CS',segType)) %>%
            dplyr::mutate(next_segType= ifelse(segName %in% class.switch.genes[-8],
                                    'CS',next_segType)) %>%
            dplyr::mutate(prev_segType = ifelse(segName %in% c(class.switch.genes,'IGHJ6'),
                                         'CS',next_segType))

          # Class switch segment guide:
          # Segments IGHA2 to IGHG3 monotonically decreasing
          # IGHG3 segment can not be lower than minimal fraction
          # (e.g. class switched B cells can not be greater than all B cells)
          # IGHM segment back to 0
          # VDJ segments then same as before


        }


        v_seg_loc <- which(segment.ranges$next_segType == 'V')
        new_seg_name <- paste0(segment.ranges$segName[v_seg_loc[1]],'_',
                               segment.ranges$segName[v_seg_loc[1] + 1])
        segment.ranges$segName[v_seg_loc[1]] <- new_seg_name
        segment.ranges <- segment.ranges[-c(v_seg_loc[1] + 1), ]
        segment.ranges <- segment.ranges %>%
          dplyr::select(segName, start = start2, end = end2) %>%
          dplyr::mutate(segName = as.character(segName))
      }else{
        segment.ranges <- vdj.segments.list[[vdj.gene]] %>%
          dplyr::select(segName = hgnc_symbol, start = start_position_hg38,
                        end = end_position_hg38) %>%
          dplyr::mutate(segName = gsub('-','_', segName)) %>%
          dplyr::filter(!is.na(start)) %>%
          dplyr::filter(!segName %in% exclude.segs) %>%
          dplyr::mutate(segType = ifelse(grepl('V',segName),'V','J')) %>%
          dplyr::mutate(next_segType = ifelse(grepl('V',dplyr::lead(segName)),'V','J')) %>%
          dplyr::mutate(prev_segType = ifelse(grepl('V',dplyr::lag(segName)),'V','J')) %>%
          dplyr::mutate(start2 = ifelse(segType == 'V',end + 1,
                                        ifelse(prev_segType == 'V',
                                               dplyr::lag(end) + 1,dplyr::lag(start)))) %>%
          dplyr::mutate(end2 = ifelse(segType == 'V',
                                      ifelse(next_segType == 'J',
                                             dplyr::lead(start) - 1,dplyr::lead(end)),
                                      start -1))

        J_seg_loc <- which(segment.ranges$next_segType == 'J')
        next_seg_name <- paste0(segment.ranges$segName[J_seg_loc[1]],'_',
                                segment.ranges$segName[J_seg_loc[1] + 1])
        segment.ranges$segName[J_seg_loc[1]] <- next_seg_name
        segment.ranges <- segment.ranges[-c(J_seg_loc[1] + 1), ]
        segment.ranges <- segment.ranges %>%
          dplyr::select(segName, start = start2, end = end2) %>%
          dplyr::mutate(segName = as.character(segName))
      }
    }else{
      if(hg19_or_38 == 'hg19'){
        if(vdj.gene %in% c('TCRG','IGH')){
          segment.ranges <- vdj.segments.list[[vdj.gene]] %>%
            dplyr::select(segName = hgnc_symbol, start = start_position_hg19,
                          end = end_position_hg19) %>%
            dplyr::mutate(segName = gsub('-','_', segName)) %>%
            dplyr::filter(!is.na(start)) %>%
            dplyr::filter(!segName %in% exclude.segs) %>%
            dplyr::mutate(segType = ifelse(grepl('V',segName),'V','J')) %>%
            dplyr::arrange(start) %>%
            dplyr::mutate(next_segType = ifelse(grepl('V',dplyr::lead(segName)),'V','J')) %>%
            dplyr::mutate(prev_segType = ifelse(grepl('V',dplyr::lag(segName)),'V','J')) %>%
            dplyr::mutate(start2 = ifelse(segType == 'J',end + 1,
                                          ifelse(prev_segType == 'J',
                                                 dplyr::lag(end) + 1,dplyr::lag(start)))) %>%
            dplyr::mutate(end2 = ifelse(segType == 'J',
                                        ifelse(next_segType == 'V',
                                               dplyr::lead(start) - 1,dplyr::lead(end)),
                                        start -1))

          V_seg_loc <- which(segment.ranges$next_segType == 'V')
          next_seg_name <- paste0(segment.ranges$segName[V_seg_loc[1]],'_',
                                  segment.ranges$segName[V_seg_loc[1] + 1])

          segment.ranges$segName[V_seg_loc[1]] <- next_seg_name
          segment.ranges <- segment.ranges[-c(V_seg_loc[1] + 1), ]
          segment.ranges <- segment.ranges %>%
            dplyr::select(segName, start = start2, end = end2) %>%
            dplyr::mutate(segName = as.character(segName))
        }else{
          segment.ranges <- vdj.segments.list[[vdj.gene]] %>%
            dplyr::select(segName = hgnc_symbol, start = start_position_hg19,
                          end = end_position_hg19) %>%
            dplyr::mutate(segName = gsub('-','_', segName)) %>%
            dplyr::filter(!is.na(start)) %>%
            dplyr::filter(!segName %in% exclude.segs) %>%
            dplyr::mutate(segType = ifelse(grepl('V',segName),'V','J')) %>%
            dplyr::arrange(start) %>%
            dplyr::mutate(next_segType = ifelse(grepl('V',dplyr::lead(segName)),'V','J')) %>%
            dplyr::mutate(prev_segType = ifelse(grepl('V',dplyr::lag(segName)),'V','J')) %>%
            dplyr::mutate(start2 = ifelse(segType == 'V',end + 1,
                                          ifelse(prev_segType == 'V',
                                                 dplyr::lag(end) + 1,dplyr::lag(start)))) %>%
            dplyr::mutate(end2 = ifelse(segType == 'V',
                                        ifelse(next_segType == 'J',
                                               dplyr::lead(start) - 1,dplyr::lead(end)),
                                        start -1))

          J_seg_loc <- which(segment.ranges$next_segType == 'J')
          next_seg_name <- paste0(segment.ranges$segName[J_seg_loc[1]],'_',
                                  segment.ranges$segName[J_seg_loc[1] + 1])
          segment.ranges$segName[J_seg_loc[1]] <- next_seg_name
          segment.ranges <- segment.ranges[-c(J_seg_loc[1] + 1), ]
          segment.ranges <- segment.ranges %>%
            dplyr::select(segName, start = start2, end = end2) %>%
            dplyr::mutate(segName = as.character(segName))
        }

      }
    }
  }



  if(GC.correct){
    if(is.null(customNorm)){
      test.logR2 <-   test.logR %>%
        dplyr::select(pos, Ratio, exon.gc, exon.gc2,smooth.gc, smooth.gc2)
    }else{
      col.loc1 <- which(colnames(test.logR) %in% c('pos','Ratio','exon.gc',
                                                   'exon.gc2','smooth.gc',
                                                   'smooth.gc2',norm_col_names))
      test.logR2 <-   test.logR[,col.loc1]
    }

  }else{
    test.logR2 <-   test.logR %>%
      dplyr::select(pos, Ratio)
  }


  # Make sure all segments have values
  segment.positions.lengths <- sapply(seq_len(dim(segment.ranges)[1]), function(x){
    start_loc <- segment.ranges$start[x]
    stop_loc <- segment.ranges$end[x]
    pos_dim <- test.logR2 %>%
      dplyr::filter(pos >= start_loc & pos <= stop_loc) %>%
      dim()
    return(pos_dim[1])
  })

  # Update the segment.ranges to remove the ones with 0 positions
  segment.ranges.upd <- segment.ranges
  if(length(which(segment.positions.lengths == 0)) > 0){
    for(i in rev(which(segment.positions.lengths == 0))){
      upd.end <- segment.ranges.upd$end[i]
      segment.ranges.upd <- segment.ranges.upd[-i, ]
      segment.ranges.upd$end[i - 1] <- upd.end
    }
  }



  # Add columns to the matrix
  for(i in seq_len(dim(segment.ranges.upd)[1])){
    test.logR2 <- create_col_seg(test.logR2,
                                 start_pos = segment.ranges.upd$start[i],
                                 end_pos = segment.ranges.upd$end[i],
                                 col.name = segment.ranges.upd$segName[i])
  }

  colnames(test.logR2) <- gsub('-','_',colnames(test.logR2))
  segment.ranges.upd$segName <- gsub('-','_',segment.ranges.upd$segName)

  # Remove IGHM column if there
  segment.ranges.names <- segment.ranges.upd$segName

  if(vdj.gene == 'IGH' & classSwitch){
    ighm.col <- which(colnames(test.logR2) == 'IGHM')
    if(length(ighm.col) > 0){
      test.logR2 <- test.logR2[,-ighm.col]
    }
    segment.ranges.names <- setdiff(segment.ranges.names, 'IGHM')
  }

  if(GC.correct){
    if(is.null(customNorm)){
      model.lm <- paste0('Ratio ~ 0 + exon.gc + exon.gc2 + smooth.gc + smooth.gc2 + ',
                         (paste0(segment.ranges.names, collapse = ' + ')))
    }else{
      model.lm <- paste0('Ratio ~ 0 + exon.gc + exon.gc2 + smooth.gc + smooth.gc2 + ',
                         paste0(norm_col_names, collapse = '+'),' + ',
                         (paste0(segment.ranges.names, collapse = ' + ')))
    }
  }else{
    model.lm <- paste0('Ratio ~ 0 + ',
                       (paste0(segment.ranges.names, collapse = ' + ')))
  }

  test.lm.unconstrained <- lm(model.lm, test.logR2)

  if(vdj.gene == 'IGH' & classSwitch){
    start_segment <- setdiff(segment.ranges.upd$segName, class.switch.genes)[1]
  }else{
    start_segment <- segment.ranges.upd$segName[1]
  }

  if(grepl('J',start_segment)){
    myConstraints1 <- segment.ranges.upd %>%
      dplyr::mutate(segName2 = paste0(dplyr::lag(segName,default = '0'),'\n')) %>%
      dplyr::mutate(comparison = ifelse(grepl('J',segName), ' < ',' > ')) %>%
      dplyr::mutate(constraint = paste0(segName, comparison, segName2))    %>%
      dplyr::select(constraint) %>% `[[`(1)
  }else{
    myConstraints1 <- segment.ranges.upd %>%
      dplyr::mutate(segName2 = paste0(dplyr::lag(segName,default = '0'),'\n')) %>%
      dplyr::mutate(comparison = ifelse(grepl('J',segName),
                                 ifelse(grepl('V',segName),' < ',' > '),' < ')) %>%
      dplyr::mutate(constraint = paste0(segName, comparison, segName2))    %>%
      dplyr::select(constraint) %>% `[[`(1)

  }

  if(vdj.gene == 'IGH' & classSwitch){
    # Fix Class switching comparisons (remove IGHM):
    # Class switch constraint loc:
    class.switch.constraint.loc <- lapply(setdiff(class.switch.genes,'IGHM'),
                                          function(x) grep(x,myConstraints1)) %>% unlist() %>% unique()
    IGHM.constraint.loc <- grep('IGHM',myConstraints1)
    class.switch.constraint.loc <- setdiff(class.switch.constraint.loc, IGHM.constraint.loc)
    if(length(class.switch.constraint.loc) == 0 | length(IGHM.constraint.loc) != 2){
      warning('Not enough bases with coverage in IGH class switch region')
      return(NULL)
    }

    myConstraints1[class.switch.constraint.loc] <- gsub('>','<',myConstraints1[class.switch.constraint.loc])

    final_cs_seg <- segment.ranges.upd$segName[max(which(segment.ranges.upd$segName %in% setdiff(class.switch.genes, 'IGHM')))]
    first_v_seg <- segment.ranges.upd$segName[min(grep('J',segment.ranges.upd$segName))]
    JV_seg <- segment.ranges.upd$segName[min(grep('V',segment.ranges.upd$segName))]

    myConstraints1[IGHM.constraint.loc[1]] <- paste0(final_cs_seg, ' > ', JV_seg, '\n')
    myConstraints1[IGHM.constraint.loc[2]] <-  paste0(first_v_seg, ' < 0\n')
  }

  myConstraints1 <- myConstraints1 %>% paste0(collapse = ' ')

  myConstraints <- paste0('\n ',myConstraints1, ' ',
                          segment.ranges.upd$segName[length(segment.ranges.upd$segName)],
                          ' < 0')

  restr.lm <- restriktor::conLM.lm(test.lm.unconstrained,
                                   constraints = myConstraints,
                                   se = 'none',mix.weights = 'none',
                                   control = list(absval = restriktor.absval))

  num.par <- length(restr.lm$parTable$est)
  num.obs <- dim(test.logR2)[1]

  vdj_segment_output <- data.frame(logR = restr.lm$b.restr) %>%
    tibble::rownames_to_column('segment') %>%
    dplyr::mutate(VDJ.tcell.fraction = 1- 2^logR) %>%
    dplyr::mutate(sample = sample_name) %>%
    dplyr::select(sample, segment, logR, VDJ.tcell.fraction) %>%
    dplyr::mutate(modelR2 =  restr.lm$R2.reduced) %>%
    dplyr::mutate(model.loglik =  restr.lm$loglik) %>%
    dplyr::mutate(model.s2=  restr.lm$s2) %>%
    dplyr::mutate(model.BIC = -2*restr.lm$loglik + log(num.obs)*num.par)

  colnames(vdj_segment_output)[4] <- ifelse(vdj.gene %in% c('TCRA','TCRB','TCRG'),
                                            paste0(vdj.gene,'.tcell.fraction'),
                                            paste0(vdj.gene,'.bcell.fraction'))

  # Remove some of the large objects
  rm(restr.lm)
  gc()

  return(vdj_segment_output)
}
