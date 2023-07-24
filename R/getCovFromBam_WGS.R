#' Function to create TCRA coverage file from bams
#'
#' @param bamPath Path to bam file
#' @param outPath Path to output directory
#' @param vdj.gene VDJ gene to use
#' @param hg19_or_38 Human genome version
#' @param samtools.loc Path to working version of samtools
#' @param index.loc Path to BAM index if in different path to BAM
#' @param customSeg Custom segment file
#' @name getCovFromBam_WGS
#' @export

getCovFromBam_WGS <- function(bamPath, outPath, vdj.gene, hg19_or_38,
                              samtools.loc = '', index.loc = NULL, customSeg = NULL){
  # Requires samtools to be installed and working!
  if(!(file.exists(bamPath))){
    stop('Can not find bam file')
  }

  vdj.chr.df <- data.frame(gene = c('TCRA','TCRB','TCRG','IGH','IGL','IGK','TCRD'),
                           chr = c('chr14','chr7','chr7','chr14','chr22','chr2','chr14'))

  if(is.null(customSeg)){
    # data("vdj_seg_list")
    seg.name <- paste0(vdj.gene, '_', hg19_or_38)
    vdj.seg <- vdj_seg_list[[seg.name]]
  }else{
    vdj.seg <- customSeg
  }


  vdj.chr <- as.character(vdj.chr.df$chr[which(vdj.chr.df$gene == vdj.gene)])

  vdj.start <- vdj.seg[vdj.seg$segName == 'all','start']
  vdj.end <- vdj.seg[vdj.seg$segName == 'all','end']

  cov.name <- gsub('.bam','',basename(bamPath))
  cov.output.files <- paste0(outPath,cov.name, '_',vdj.gene,'.txt')

  # Check does bai file exist
  if(!is.null(index.loc)){
    bai.path <- index.loc
  }else{
    bai.path <- paste0(bamPath,'.bai')
    bai.path2 <- paste0(gsub('.bam','',bamPath),'.bai')
    crai.path <- paste0(bamPath,'.crai')
    crai.path2 <- paste0(gsub('.cram','',bamPath),'.crai')
    if(!(file.exists(bai.path) | file.exists(bai.path2) | file.exists(crai.path) | file.exists(crai.path2))){
      stop('No index bai or crai file found for bam or cram, please index first before preceding')
    }
  }


  # Check if bam has chr or not before
  samtools.chr.check <- paste0(samtools.loc,'samtools idxstats ', bamPath,' | head -n 2')
  chr.check.output <- system(samtools.chr.check, intern = TRUE, ignore.stderr = TRUE)[2]
  chr.check.output2 <- strsplit(chr.check.output,'\t')[[1]][1]
  chr.present <- grepl('chr',chr.check.output2)

  chr.to.use <- ifelse(chr.present, vdj.chr ,gsub('chr','',vdj.chr))

  if(!is.null(index.loc)){
    samtools.cmds <- paste0(samtools.loc,'samtools depth ',bamPath,
                            ' -X ', bai.path,
                            ' -q 20 -Q 20 -r ',chr.to.use,':',vdj.start,'-',vdj.end,' > ',
                            cov.output.files)
  }else{
    samtools.cmds <- paste0(samtools.loc,'samtools depth ',bamPath,
                            ' -q 20 -Q 20 -r ',chr.to.use,':',vdj.start,'-',vdj.end,' > ',
                            cov.output.files)
  }


  # Replace with processx!
  sapply(samtools.cmds, system)

  return(cov.output.files)
}
