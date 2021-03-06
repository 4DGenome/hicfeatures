% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{hg38_NcoI_10000}
\alias{hg38_NcoI_10000}
\title{Genomic features for hg38 genome and NcoI restriction enzyme at 10 Kbp}
\format{A data frame with 303113 rows and 5 variables:
 \describe{
   \item{chr:}{chromosome}
   \item{map:}{mappability as computed by gem}
   \item{res:}{restriction enzyme density per 1 Kbp computed by Biostrings::matchPattern()}
   \item{cg:}{cg content as computed by bedtools}
   \item{bin:}{genomic bin with the format chromosome:start_position}
   \item{pos:}{start postion of the genomic bin}
 }}
\usage{
hg38_NcoI_10000
}
\description{
A \code{data.frame} containing the mappability, restriction
enzyme density and CG proportion of the hg38 genome and
NcoI restriction enzyme in 10 Kbp bins
}
\keyword{datasets}
