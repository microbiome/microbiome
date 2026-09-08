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

    # Compute the abundance matrix once and take the margins from it, so
    # that the summary works for any supported backend. This also avoids
    # recomputing the same sums a dozen times.
    a <- abundances(x)
    ssums <- colSums(a)
    tsums <- rowSums(a)

    ave <- sum(ssums)/ncol(a)

    comp <- length(which(ssums > 1))

    if (comp == 0)
    {
        message("Compositional = YES", fill = 2)
    } else
    {
        message("Compositional = NO", fill = 2)
    }


    minR <- paste0("1] Min. number of reads = ", min(ssums))
    maxR <- paste0("2] Max. number of reads = ", max(ssums))
    tR <- paste0("3] Total number of reads = ", sum(ssums))
    aR <- paste0("4] Average number of reads = ", ave)
    mR <- paste0("5] Median number of reads = ", median(ssums))

    if (any(tsums <= 1) == TRUE)
    {
        sR <- paste0("6] Any OTU sum to 1 or less? ", "YES")
    } else
    {
        sR <- paste0("6] Any OTU sum to 1 or less? ", "NO")

    }

    zno <- paste0("7] Sparsity = ", length(which(a == 0))/length(a))

    sR1 <- paste0("8] Number of singletons = ",
        length(tsums[tsums <= 1]))
    sR2 <- paste0("9] Percent of OTUs that are singletons
        (i.e. exactly one read detected across all samples)",
        mean(tsums == 1) * 100)
        #length(taxa_sums(x)[taxa_sums(x) <= 1])/nrow(otu_table(x)) * 
        #    100)
    svar <- paste0("10] Number of sample variables are: ", ncol(meta(x)))
    sam_var <- colnames(meta(x))

    message(minR, maxR, tR, aR, mR, zno, sR, sR1, sR2, svar, sam_var, 
        fill = 2)

    # Return?
    list(minR, maxR, tR, aR, mR, zno, sR, sR1, sR2, svar, sam_var)

}
