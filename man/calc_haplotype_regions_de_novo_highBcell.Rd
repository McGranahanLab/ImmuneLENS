% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_haplotype_regions_de_novo_highBcell.R
\name{calc_haplotype_regions_de_novo_highBcell}
\alias{calc_haplotype_regions_de_novo_highBcell}
\title{Calculate predicted IGH germline haplotype assuming possible high B cell content}
\usage{
calc_haplotype_regions_de_novo_highBcell(
  input_cov,
  kb_len_threshold = 5,
  round_solution = TRUE
)
}
\arguments{
\item{input_cov}{data frame of positions and coverage values for IGH loci}

\item{kb_len_threshold}{threshold for length of genomic regions used in normalisation method (default = 5kb)}

\item{round_solution}{default = TRUE}
}
\value{
regions used for later normalisation
}
\description{
Calculate predicted IGH germline haplotype assuming possible high B cell content
}
