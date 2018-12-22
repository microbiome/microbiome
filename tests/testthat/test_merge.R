context("Merging operations")

test_that("Merging ok", {

  data(dietswap)
  ps <- dietswap
  expect_equal(ntaxa(merge_taxa2(ps, pattern = "^Prevotella")), 127)
  expect_equal(ntaxa(merge_taxa2(ps, taxa = taxa(ps)[grep("Prevotella", taxa(ps))])), 127)  
  expect_true("Prevotella" %in% taxa(merge_taxa2(ps, pattern = "^Prevotella", name = "Prevotella")))
  expect_equal(abundances(aggregate_taxa(ps, level = "Phylum"))["Firmicutes", 2],
               sum(abundances(ps)[rownames(tax_table(ps))[which(tax_table(ps)[, "Phylum"] == "Firmicutes")],2]))


  #x <- dietswap
  library(phyloseq)
  x <- GlobalPatterns
  a1 <- abundances(aggregate_taxa(x, level = "Genus"))
  a2 <- colSums(abundances(x)[which(tax_table(x)[, "Genus"] == "Dialister"), ])
  expect_equal(a1["Dialister",], a2)

  # Top taxa aggregation
  x <- dietswap
  a3 <- abundances(aggregate_top_taxa(x, top = 3, level = "Phylum"))
  x1 <- a3["Bacteroidetes",]
  inds <- which(tax_table(x)[, "Phylum"] == "Bacteroidetes")
  mytaxa <- rownames(abundances(x))[inds]
  x2 <- colSums(abundances(x)[inds,])
  expect_equal(x1, x2)    
 

})

