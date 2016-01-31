core.sum <- function(data, intTr, prevalenceTr) {

    # Prevalence for each OTU
    prevalences <- rowSums(data > intTr)

    # Number of OTUs above a given prevalence threshold
    nOTUs <- sum(prevalences >= prevalenceTr)

    nOTUs

}


