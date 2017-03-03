context("diversity")

test_that("diversity works correctly", {

  data(peerj32)
  pseq <- peerj32$phyloseq
  
  # For phyloseq
  expect_true(is.vector(diversity(pseq, split = FALSE)))
  expect_true(is.data.frame(diversity(pseq, split = TRUE)))  
  
})

