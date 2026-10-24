test_that("camargo evenness is computed over observed taxa (#231)", {

  # Reference implementation of Camargo (1992)
  ref <- function(x) {
    p <- x[x > 0] / sum(x)
    1 - sum(abs(outer(p, p, "-"))) / 2 / length(p)
  }

  # Number of leading zeroes must not matter
  for (k in c(2, 5, 10)) {
    x <- c(rep(0, k), 0.99, rep(5e-05, 5), rep(0, 110))
    expect_equal(evenness(x, "camargo")$camargo, ref(x))
  }

  # Order invariance
  expect_equal(evenness(c(0, 0, 0, 5, 1, 1), "camargo")$camargo,
               evenness(c(5, 1, 1, 0, 0, 0), "camargo")$camargo)

  # Perfectly even community
  expect_equal(evenness(c(2, 2, 2, 0, 2), "camargo")$camargo, 1)

  # Hand-computed values
  expect_equal(evenness(c(3, 1), "camargo")$camargo, 0.75)
  expect_equal(evenness(c(9, 1), "camargo")$camargo, 0.6)

  # Single taxon and empty samples
  expect_equal(evenness(c(0, 4, 0), "camargo")$camargo, 1)
  expect_true(is.na(evenness(c(0, 0, 0), "camargo")$camargo))

  # Matrix input agrees with the reference
  data(dietswap)
  otu <- abundances(dietswap)[, 1:10]
  expect_equal(unname(evenness(otu, "camargo")$camargo),
               unname(apply(otu, 2, ref)))

})
