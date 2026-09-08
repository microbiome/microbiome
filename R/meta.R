#' @title Retrieve Sample Metadata as Data Frame
#' @description The output of the phyloseq::sample_data() function does not
#' return data.frame, which is needed for many applications.
#' This function retrieves the sample data as a data.frame. For
#' SummarizedExperiment-derived objects (including
#' TreeSummarizedExperiment) the colData is returned instead.
#' @param x a phyloseq or SummarizedExperiment-derived object
#' @return Sample metadata as a data.frame
#' @export
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples data(dietswap); df <- meta(dietswap)
#' @seealso \code{\link{sample_data}} in the \pkg{phyloseq} package
#' @keywords utilities
meta <- function(x) {

    if (.is_se(x)) {

        cd <- SummarizedExperiment::colData(x)

        # A colData with no columns coerces to a zero-row data.frame,
        # which would then reject the rownames. Build it explicitly so
        # that the sample count is preserved either way.
        if (ncol(cd) == 0) {
            return(data.frame(row.names=colnames(x)))
        }

        df <- as.data.frame(cd)
        rownames(df) <- colnames(x)

        return(df)

    }

    .deprecate_phyloseq(x)

    df <- as(sample_data(x), "data.frame")
    rownames(df) <- sample_names(x)
    df
}


