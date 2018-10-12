#' Prepare genomic features per genomic bin
#'
#' This function takes a genome, enzyme and resolution and produces a \code{data.frame} with the corresponding genomic features. It's just a wrapper to the \code{data()} contained in the package. In case that the desired resolution is not available in the pre-computed data, it interpolates the genomic features from the nearest larger resolution available.
#' @import dplyr
#' @param geno A \code{string} encoding a genome assemply. Currently, only \code{hg38}, \code{mm10},  \code{dm3} and \code{dm6} are supported.
#' @param geno A \code{string} encoding a restriction enzyme. Currently, only \code{BspHI}, \code{DpnII}, \code{HindIII}, \code{HinfI}, \code{MboI}, \code{MspI}, \code{NcoI} and \code{NlaIII} are supported.
#' @param reso A \code{integer} with the desired bin size (a.k.a. resolution).
#' @return A \code{data.frame} containing 6 variables:
#'  \describe{
#'    \item{chr:}{chromosome}
#'    \item{map:}{mappability as computed by gem}
#'    \item{res:}{restriction enzyme density per 1 Kbp computed by Biostrings::matchPattern()}
#'    \item{cg:}{cg content as computed by bedtools}
#'    \item{bin:}{genomic bin with the format chromosome:start_position}
#'    \item{pos:}{start postion of the genomic bin}
#'  }
#' 
#' @export
#' @examples
#' plot(0)

get_genomic_features <- function(geno, enz, reso){

    enz <- ifelse(enz == "MboI", "DpnII", enz)
    av_genomes <- c("dm3", "dm6", "hg38", "mm10")
    av_enzymes <- c("BspHI", "DpnII", "HindIII", "HinfI", "MspI", "NcoI", "NlaIII")
    av_resolutions <- c(1000000, 100000, 10000, 1000, 5000000, 500000, 50000, 5000)
    
    stopifnot(geno %in% av_genomes)
    stopifnot(enz %in% av_enzymes)
    
    if(reso %in% av_resolutions){
        out <- paste0(geno,
                      "_", enz, "_",
                      as.integer(reso)) %>%
            get()
        return(out)
    }

    aux <- (reso - av_resolutions)
    aux <- av_resolutions[aux >= 0]
    best_reso <- max(aux)

    out <- paste0(geno,
                  "_", enz, "_",
                  as.integer(best_reso)) %>%
        get() %>%
        mutate(pos = floor(pos /reso) * reso,
               pos = as.integer(pos),
               bin = paste0(chr, ":", pos)) %>%
        group_by(chr, pos, bin) %>%
        summarize(res = mean(res, na.rm = T),
                  cg = mean(cg, na.rm = T),
                  map = mean(map, na.rm = T)) %>%
        ungroup() %>%
        dplyr::select(chr, map, res, cg, bin, pos)

    out
    
}
