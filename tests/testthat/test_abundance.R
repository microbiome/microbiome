

context("abundance")

test_that("abundance picking works correctly", {

  data(dietswap)
  pseq <- dietswap
  
  expect_true(is.matrix(abundances(pseq)))
  expect_true(is.vector(abundances(pseq)[,1]))  
  expect_true(is.matrix(abundances(pseq, transform = "compositional")))  


})

