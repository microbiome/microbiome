# For more details on unit testing, see
# http://master.bioconductor.org/developers/how-to/unitTesting-guidelines/

context("bimodality")

test_that("bimodality scores work correctly", {

  # Two equally frequent values give 1, a uniform distribution 5/9
  expect_equal(bimodality_sarle(rep(c(0, 1), 50), type = "Sarle.asymptotic"), 1)
  expect_equal(bimodality_sarle(seq(0, 1, length.out = 1e5),
      type = "Sarle.asymptotic"), 5/9, tolerance = 1e-4)
  expect_equal(bimodality_sarle(seq(100), type = "Sarle.finite.sample"), 0.5283013, tolerance=1.0e-6)
  
})

