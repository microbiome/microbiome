context("compositionality")

test_that("Compositionality tests are correct", {

  data(dietswap)
  pseq <- dietswap
  expect_true(is.compositional(transform(pseq, "compositional")))
  expect_false(is.compositional(pseq))
  
})

