#' @title Summarize phyloseq object
#' @description Prints basic information of data.
#' @details The summarize_phyloseq function will give information on weather
#'    data is compositional or not, reads (min. max, median, average),
#'    sparsity, 
#' presence of singletons and sample variables. 
#' @param x Input is a \code{\link{phyloseq-class}} object.
#' @return Prints basic information of \code{\link{phyloseq-class}} object.
#' @export
#' @examples
#' data(dietswap)
#' summarize_phyloseq(dietswap)
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@@gmail.com}
#' @keywords utilities
#' 
summarize_phyloseq <- function(x)
{
    
    ave <- minR <- maxR <- tR <- aR <- mR <- sR <- sR1 <- sR2 <- svar <- NULL
    sam_var <- zno <- comp <- NULL
    
    ave <- sum(sample_sums(x))/nsamples(x)
    
    comp <- length(which(colSums(abundances(x)) > 1))
    
    if (comp == 0)
    {
        message("Compositional = YES", fill = 2)
    } else
    {
        message("Compositional = NO", fill = 2)
    }
    
    
    minR <- paste0("1] Min. number of reads = ", min(sample_sums(x)))
    maxR <- paste0("2] Max. number of reads = ", max(sample_sums(x)))
    tR <- paste0("3] Total number of reads = ", sum(sample_sums(x)))
    aR <- paste0("4] Average number of reads = ", ave)
    mR <- paste0("5] Median number of reads = ", median(sample_sums(x)))
    
    if (any(taxa_sums(x) <= 1) == TRUE)
    {
        sR <- paste0("6] Any OTU sum to 1 or less? ", "YES")
    } else
    {
        sR <- paste0("6] Any OTU sum to 1 or less? ", "NO")
        
    }
    
    zno <- paste0("7] Sparsity = ", length(which(abundances(x) == 
        0))/length(abundances(x)))
    
    sR1 <- paste0("8] Number of singletons = ",
        length(taxa_sums(x)[taxa_sums(x) <= 
        1]))
    sR2 <- paste0("9] Percent of OTUs that are singletons 
        (i.e. exactly one read detected across all samples)",
        mean(taxa_sums(x) == 1) * 100)
        #length(taxa_sums(x)[taxa_sums(x) <= 1])/nrow(otu_table(x)) * 
        #    100)
    svar <- paste0("10] Number of sample variables are: ", ncol(meta(x)))
    sam_var <- colnames(meta(x))
    
    message(minR, maxR, tR, aR, mR, zno, sR, sR1, sR2, svar, sam_var, 
        fill = 2)

    # Return?
    list(minR, maxR, tR, aR, mR, zno, sR, sR1, sR2, svar, sam_var)
    
}
