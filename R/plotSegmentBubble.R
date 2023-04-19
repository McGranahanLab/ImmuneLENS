#' Function to visualise segment output as bubble plot
#'
#' @param segment_out data frame output from runTcellExtrect_WGS containing segment fractions
#' @param V_or_J plot V or J segment proportions
#' @param text_size size of text to plot in bubbles
#' @return Bubble plots
#' @importFrom ggplot2 ggplot aes geom_polygon geom_text coord_equal theme_bw theme element_blank
#' @name plotSegmentBubble
#' @export

plotSegmentBubble <- function(segment_out, V_or_J, text_size =  6){
  if(!V_or_J %in% c('V','J')){
    stop('V_or_J must be either "V" or "J"')
  }
  segment.fraction <- segment.prop <- x <- y <- id <- v_seg <- segment <- NULL

  bubble.df <- segment_out %>%
    dplyr::filter(grepl(V_or_J, segment)) %>%
    dplyr::mutate(segment.prop = (segment.fraction)/sum(segment.fraction)) %>%
    dplyr::filter(segment.prop > 1e-4)

  limits = c(-1,1)
  res <- packcircles::circleRepelLayout(bubble.df$segment.prop,
                                        xlim = limits, ylim = limits)
  res$layout$v_seg <- bubble.df$segment
  dat <- packcircles::circlePlotData(res$layout,id.col = 'v_seg')

  p1 <- ggplot(data = dat, aes(x, y)) +
    geom_polygon(aes(fill = id)) +
    geom_text(aes(label = v_seg),data = res$layout, size = text_size/2.845276) +
    coord_equal(xlim=c(min(dat$x), max(dat$x)), ylim=c(min(dat$y),max(dat$y))) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_blank(),
          legend.position = 'none',
          panel.border = element_blank())

  print(p1)

}





