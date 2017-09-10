#' @title RDA for phyloseq objects
#' @description RDA for phyloseq objects based on the \code{\link{rda}}
#'    function from the \pkg{vegan} package.
#' @param x \code{\link{phyloseq-class}} object
#' @param y Variable to apply in RDA visualization.
#' @param scale See help(rda)
#' @param na.action See help(rda)
#' @param ... Other arguments to be passed
#' @return rda result. See help(vegan::rda)
#' @export
#' @examples
#'   data(peerj32) # Data from https://peerj.com/articles/32/
#'   pseq <- peerj32$phyloseq
#'   pseq.trans <- transform(pseq, "hell") # Hellinger transform
#'   # rda.result <- rda_pseq(pseq.trans, "time", scale = TRUE)
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
rda_pseq <- function (x, y, scale = FALSE, na.action = na.fail, ...) {

  # Microbiota profiling data (44 samples x 130 bacteria)
  otu <- t(abundances(x))

  # Pick the indicated annotation field
  if (!y %in% sample_variables(x)) {    
      stop(paste("The variable y ('", y, "') is not available in the phyloseq object i.e. sample_data(x). Only use variables listed in sample_variables(x) ie. o
ne of the following: ", paste(names(sample_data(x)), collapse = " / "), sep = ""))
   }

  if (!"sample" %in% sample_variables(x)) {
    sample_data(x)$sample <- rownames(sample_data(x))
  }

  annot <- factor(sample_data(x)[[y]])
  names(annot) <- sample_data(x)$sample

  r <- vegan::rda(otu ~ annot, scale = scale, na.action = na.action) 
  
  r

}



