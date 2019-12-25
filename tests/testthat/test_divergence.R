context("divergence")

test_that("Divergence", {


  data(dietswap)
  pseq <- subset_samples(dietswap, nationality == 'AFR')
  reference <- apply(abundances(pseq), 1, median)

  expect_equal(divergence(reference, reference, method = "bray"), 0)
  expect_equal(divergence(abundances(pseq)[,1], reference, method = "bray"), 0.3449469, tolerance = 1e-6)  

  b <- divergence(pseq, reference, method = "bray")
   
})

