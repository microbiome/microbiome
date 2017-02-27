#' @title Retrieve Phyloseq Metadata as Data Frame
#' @description The output of the phyloseq::sample_data() function does not return data.frame, which is needed for many applications. This function retrieves the sample data as a data.frame
#' @param x a phyloseq object
#' @export
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples data(dietswap); df <- meta(dietswap)
#' @keywords utilities
pseq_metadata <- function (x) {

  .Deprecated(new = meta)
  meta(x)

}
