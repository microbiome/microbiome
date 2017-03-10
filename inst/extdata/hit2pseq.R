#' @title HITChip to phyloseq
#' @description Convert HITChip data into phyloseq format.
#' @param otu Sample x OTU absolute HITChip signal
#' @param meta Sample x features metadata data.frame
#' @param taxonomy OTU x Taxonomy data.frame (HITChip taxonomy used by default)
#' @param detection.limit HITChip signal detection limit (absence / presence)
#' @return phyloseq object
#' @importFrom phyloseq otu_table
#' @importFrom phyloseq tax_table
#' @importFrom phyloseq phyloseq
#' @importFrom phyloseq sample_data
#' @importFrom phyloseq merge_phyloseq
#' @examples 
#'  \dontrun{
#'   library(microbiome)
#'   data(peerj32)
#'   otu <- peerj32$microbes
#'   meta <- peerj32$meta
#'   pseq <- hitchip2physeq(otu, meta)
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

  # Create phyloseq object
  pseq <- phyloseq(OTU)

  # --------------------------

  # Construct taxonomy table
  # If nrow(otumat) then it is probe-level data and no taxonomy should be given
  # for that by default
  if (is.null(taxonomy) && nrow(otumat) < 3000) {

    #ph <- GetPhylogeny("HITChip")
    ph = get_hitchip_taxonomy("HITChip", phylogeny.version = "full", data.dir = NULL)
    ph <- unique(ph[, c("L1", "L2", "species")])
    colnames(ph) <- c("Phylum", "Genus", "Phylotype")
    taxonomy <- ph
    input.level <- colnames(ph)[[which.max(apply(ph, 2, function (x) {sum(rownames(otumat) %in% x)}))]]
    if (input.level == "Genus") {
      taxonomy <- unique(taxonomy[, c("Phylum", "Genus")])
    } else if (input.level == "Phylum") {
      taxonomy <- data.frame(Phylum = unique(taxonomy[, c("Phylum")]))
    }
    rownames(taxonomy) <- as.character(taxonomy[[input.level]])
  }

  if (!all(rownames(otumat) %in% rownames(taxonomy)) && nrow(otumat) < 3000) {
      warning(paste("Some OTUs are missing from the taxonomy tree!", paste(setdiff(rownames(otumat), rownames(taxonomy)), collapse = " / ")))
      # Common probes or OTUs
      coms <- intersect(rownames(otumat), rownames(taxonomy))
      # Only keep probes that have taxonomy information
      otumat <- otumat[coms, ]
      taxonomy <- taxonomy[coms, ]
   }

   if (!is.null(taxonomy) || nrow(otumat) < 3000) {
     TAX <- tax_table(as.matrix(taxonomy[rownames(otumat), ]))

     # Combine OTU and Taxon matrix into Phyloseq object
     pseq <- merge_phyloseq(pseq, TAX)
  }
  
  # -------------------------

  if (!is.null(meta)) {
  
    # Metadata
    rownames(meta) <- as.character(meta$sample)
    sampledata <- sample_data(meta[colnames(otumat),])
    pseq <- merge_phyloseq(pseq, sampledata)
    
  }
  
  # --------------------------

  # We could also add phylotree between OTUs
  # source("tree.R")
  # pseq <- merge_phyloseq(pseq, tree2)
  # pseq <- merge_phyloseq(pseq, random_tree)

  # --------------------------

  pseq
 
}


