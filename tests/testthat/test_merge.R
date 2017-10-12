context("Merging operations")

test_that("Merging ok", {

  data(dietswap)
  ps <- dietswap
  expect_equal(ntaxa(merge_taxa2(ps, pattern = "^Prevotella")), 127)
  expect_equal(ntaxa(merge_taxa2(ps, taxa = taxa(ps)[grep("Prevotella", taxa(ps))])), 127)  
  expect_true("Prevotella" %in% taxa(merge_taxa2(ps, pattern = "^Prevotella", name = "Prevotella")))

})

