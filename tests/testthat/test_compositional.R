context("compositionality")

test_that("Compositionality tests are correct", {

  data(dietswap)
  pseq <- dietswap
  expect_false(is_compositional(pseq))  
  expect_true(is_compositional(transform(pseq, "compositional")))

  
})

