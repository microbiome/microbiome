core.sum <- function(data, intTr, prevalenceTr) {

    d.bin <- data > intTr
    prevalences <- rowSums(d.bin)
    nOTUs <- sum(prevalences >= prevalenceTr)
    return(nOTUs)

}


