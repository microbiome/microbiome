#' hitchip2physeq
#'
#' Convert HITChip data into phyloseq format
#'
#' @param data Sample x OTU absolute HITChip signal
#' @param meta Sample x features metadata data.frame
#' @param taxonomy OTU x Taxonomy data.frame (HITChip taxonomy used by default)
#' @return phyloseq object
#' @import phyloseq
#'
#' @examples 
#'   library(microbiome)
#'   data(peerj32)
#'   data <- peerj32$microbes
#'   meta <- peerj32$meta
#'   physeq <- hitchip2physeq(data, meta)
#'
#' @export
#' @references Utilizes the phyloseq package, see citation("phyloseq"). 
#'             For this function, see citation('microbiome').  
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
hitchip2physeq <- function (data, meta, taxonomy = NULL) {

  # OTU matrix	       
  # OTU x Sample: absolute 'read counts'
  x <- t(data) - 10^1.8 # Detection limit
  x[x < 0] <- 0
  x <- 1 + x

  # Discretize to get 'counts'
  otumat <- round(x)
  rownames(otumat) <- gsub("Clostridium \\(sensu stricto\\)", "Clostridiales", rownames(otumat))
  OTU <- otu_table(otumat, taxa_are_rows = TRUE)

  # --------------------------

  # Construct taxonomy table
  if (is.null(taxonomy)) {
    # Assuming for now that the input data is L2 level
    ph <- GetPhylogeny("HITChip")
    ph <- unique(ph[, c("L1", "L2")])
    ph$L2 <- gsub("Clostridium \\(sensu stricto\\)", "Clostridiales", ph$L2)
    ph <- unique(ph[, c("L1", "L2")])
    colnames(ph) <- c("Phylum", "Genus")
    taxonomy <- ph
    rownames(taxonomy) <- as.character(taxonomy$Genus)
  }

  if (!all(rownames(otumat) %in% rownames(taxonomy))) {stop("Some OTUs are missing from the taxonomy tree!")}

  TAX <- tax_table(as.matrix(taxonomy[rownames(otumat), ]))

  #----------------------------------------------

  # Combine OTU and Taxon matrix into Phyloseq object
  physeq <- phyloseq(OTU, TAX)

  # --------------------------

  # Metadata
  sampledata <- sample_data(meta[colnames(otumat),])
  physeq <- merge_phyloseq(physeq, sampledata)

  # --------------------------

  # We could also add phylotree between OTUs
  # source("tree.R")
  # physeq <- merge_phyloseq(physeq, tree2)
  # physeq <- merge_phyloseq(physeq, random_tree)

  # --------------------------

  physeq
 
}

