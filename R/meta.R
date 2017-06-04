#' @title Retrieve Phyloseq Metadata as Data Frame
#' @description The output of the phyloseq::sample_data() function does not
#' return data.frame, which is needed for many applications.
#' This function retrieves the sample data as a data.frame
#' @param x a phyloseq object
#' @return Sample metadata as a data.frame
#' @export
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples data(dietswap); df <- meta(dietswap)
#' @seealso \code{\link{sample_data}} in the \pkg{phyloseq} package
#' @keywords utilities
meta <- function(x) {
    df <- as(sample_data(x), "data.frame")
    rownames(df) <- sample_names(x)
    df
}


