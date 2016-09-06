#' @title Pick Baseline Timepoint Samples
#' @description Identify and select the baseline timepoint samples in a \class{\link{phyloseq}} object.
#' @param x \class{\link{phyloseq}} object. Assuming that the sample_data(x) has the fields "time", "sample" and "subject"
#' @param na.omit Logical. Ignore samples with no time point information. If this is FALSE, the first sample for each subject is selected even when there is no time information.
#' @return Phyloseq object with only baseline time point samples selected.
#' @details Arranges the samples by time and picks the first sample for each subject.
#' @examples
#'  data(atlas1006)
#'  a <- pick_baseline(atlas1006)
#' @export
#' @references See citation("microbiome")
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
pick_baseline <- function (x, na.omit = TRUE) {

  # Arrange samples by time, and pick the first sample for each subject
  m <- sample_data(x)

  if (!all(c("time", "sample", "subject") %in% names(m))) {
    stop("The phyloseq sample_data(x) should contain the following fields: time, sample, subject.")
  }

  m <- m %>% arrange(time)

  if (na.omit) {
    m <- m %>% filter(!is.na(time))
  }

  s <- m$sample[match(unique(m$subject), m$subject)]
  x <- subset_samples(x, sample %in% s)

  x
}
