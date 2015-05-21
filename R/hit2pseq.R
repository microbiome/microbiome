#' hitchip2physeq
#'
#' Convert HITChip data into phyloseq format
#'
#' @param otu Sample x OTU absolute HITChip signal
#' @param meta Sample x features metadata data.frame
#' @param taxonomy OTU x Taxonomy data.frame (HITChip taxonomy used by default)
#' @param detection.limit HITChip signal detection limit (absence / presence)
#' @return phyloseq object
#' @import phyloseq
#'
#' @examples 
#'  \dontrun{
#'   library(microbiome)
#'   data(peerj32)
#'   otu <- peerj32$microbes
#'   meta <- peerj32$meta
#'   physeq <- hitchip2physeq(otu, meta)
#' }
#' @export
#' @references Utilizes the phyloseq package, see citation("phyloseq"). 
#'             For this function, see citation('microbiome').  
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
hitchip2physeq <- function (otu, meta, taxonomy = NULL, detection.limit = 10^1.8) {

  # OTU x Sample matrix: absolute 'read counts'
  x <- t(otu) - detection.limit # HITChip detection limit
  x[x < 0] <- 0
  x <- 1 + x

  # -------------------------

  # Discretize to get 'counts'
  otumat <- round(x)
  OTU <- otu_table(otumat, taxa_are_rows = TRUE)

  # --------------------------

  # Construct taxonomy table
  if (is.null(taxonomy)) {

    ph <- as.data.frame(GetPhylogeny("HITChip")@.Data)
    ph <- unique(ph[, c("L1", "L2", "species")])
    colnames(ph) <- c("Phylum", "Genus", "Phylotype")
    taxonomy <- ph
    input.level <- colnames(ph)[[which.max(apply(ph, 2, function (x) {sum(rownames(otumat) %in% x)}))]]
    if (input.level == "Genus") {
      taxonomy <- unique(taxonomy[, c("Phylum", "Genus")])
    } else if (input.level == "Phylum") {
      taxonomy <- as.data.frame(Phylum = unique(taxonomy[, c("Phylum")]), ncol = 1)
      	       
    }
    rownames(taxonomy) <- as.character(taxonomy[[input.level]])

  }

  if (!all(rownames(otumat) %in% rownames(taxonomy))) {
    warning(paste("Some OTUs are missing from the taxonomy tree!", paste(setdiff(rownames(otumat), rownames(taxonomy)), collapse = " / ")))
    # Common probes
    coms <- intersect(rownames(otumat), rownames(taxonomy))
    # Only keep probes that have taxonomy information
    otumat <- otumat[coms, ]
    taxonomy <- taxonomy[coms, ]
  }

  TAX <- tax_table(as.matrix(taxonomy[rownames(otumat), ]))

  #----------------------------------------------

  # Combine OTU and Taxon matrix into Phyloseq object
  physeq <- phyloseq(OTU, TAX)

  # --------------------------

  # Metadata
  sampledata <- sample_data(meta[colnames(otumat),])
  physeq <- merge_phyloseq(physeq, sampledata)

  # Harmonize the fields
  physeq@sam_data <- harmonize_fields(physeq@sam_data)

  # --------------------------

  # We could also add phylotree between OTUs
  # source("tree.R")
  # physeq <- merge_phyloseq(physeq, tree2)
  # physeq <- merge_phyloseq(physeq, random_tree)

  # --------------------------

  physeq
 
}

