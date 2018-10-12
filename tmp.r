
# dependencies

library("dplyr")
library("stringr")
library("tidyr")

info <- list.files("data") %>%
    str_replace(".rda", "") %>%
    data.frame(info = .) %>%
    separate(info, c("genome", "enzyme", "reso"), sep = "_") %>%
    mutate(reso = as.integer(reso))

info_unique <- lapply(info, unique)

av_genomes <- info_unique$genome
av_enzymes <- info_unique$enzyme
av_resos <- info_unique$reso


enz <- "DpnII"
geno <- "hg38"
reso <- 20e3

get_genomic_features <- function(geno, enz, reso){

    enz <- ifelse(enz == "MboI", "DpnII", enz)
    av_genomes <- c("dm3", "hg38", "mm10")
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
        ungroup()

    out
    
}



get_genomic_features("hg38", "DpnII", 150000) %>%
    str
