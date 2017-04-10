# Take vector, matrix, or phyloseq and return matrix (taxa x samples)
pick_data <- function (x, compositional = FALSE) {
  # Pick the OTU data
  if (is.phyloseq(x)) {
    otu <- abundances(x)
  } else if (is.vector(x)) {
    otu <- as.matrix(x, ncol = 1)    
  } else {
    # Assuming taxa x samples count matrix
    otu <- x
  }

  if (compositional) {
    otu <- apply(otu, 2, function (x) {x/sum(x, na.rm = TRUE)})
  }

  otu
}