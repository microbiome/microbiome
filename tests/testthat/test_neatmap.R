context("neatmap sorting")

test_that("neatmap sorting works correctly", {

  data(peerj32)
  pseq <- core(subset_samples(peerj32$phyloseq, group == "Placebo"),
      detection = 5, prevalence = 50/100)
  
  # For phyloseq
  expect_true("sample-9" %in% neatsort(pseq, target = "sites"))
  expect_true("Lactococcus" %in% neatsort(pseq, target = "species"))  

  # For matrices
  expect_true("sample-9" %in% neatsort(abundances(pseq), target = "cols"))
  expect_true("Lactococcus" %in% neatsort(abundances(pseq), target = "rows"))  
  
})

