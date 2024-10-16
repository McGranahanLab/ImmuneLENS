#' Calculate median values of coverage within exons in rolling windows
#'
#' @param vdj_logR_input data frame of positions and coverage values
#' @param vdj_seg Segments used for normalisation
#' @param GCcorrect Use GC corrected output or not
#' @importFrom stats sd 
#' @return Adjusted baseline logR dataframe
#' @name baselineAdj

baselineAdj <- function(vdj_logR_input, vdj_seg, GCcorrect = TRUE){

  pos <- reads <- NULL

  ratio.col <- ifelse(GCcorrect, 'Ratio.gc.correct','Ratio')
  ratio.col <- rlang::sym(ratio.col)

  adjust.baseline.value <- vdj_logR_input %>%
    dplyr::filter((pos >= vdj_seg[3,'start'] & pos <= vdj_seg[3,'end']) |
             (pos >= vdj_seg[4,'start'] & pos <= vdj_seg[4,'end'])) %>%
    dplyr::summarise(gc.adjust = mean(!!ratio.col,na.rm = TRUE),
              CI.95.range = 1.96*sd(!!ratio.col,na.rm = TRUE)/sqrt(length(!!ratio.col)))

  # Look at focal too in GAM model

  if(GCcorrect){
    adjust.model <- mgcv::gam(Ratio.gc.correct~s(pos, bs = 'cs'), data = vdj_logR_input)
  }else{
    adjust.model <- mgcv::gam(Ratio~s(pos, bs = 'cs'), data = vdj_logR_input)
  }

  # Define the sequence of positions for prediction
  pos_sequence <- seq(vdj_seg[2, 2], vdj_seg[2, 3], by = 100)
  
  # Predict values for the sequence
  adjust.fit.values <- mgcv::predict.gam(adjust.model, newdata = data.frame(pos = pos_sequence),
                                      type = "response")
  
  # Filter the predictions to be within the desired range
  fit.loc <- which(pos_sequence > vdj_seg[2, 2] & pos_sequence < vdj_seg[2, 3])
  
  # Calculate the mean of the predictions within the range
  adjust.baseline.value2 <- list(mean(adjust.fit.values[fit.loc], na.rm = TRUE))

  adjust.value <- max(adjust.baseline.value[[1]], adjust.baseline.value2[[1]])


  vdj_logR_input[[ratio.col]] <- vdj_logR_input[[ratio.col]] - adjust.value
  return(list(vdj_logR_input))
}

