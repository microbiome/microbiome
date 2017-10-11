# For more details on unit testing, see
# http://master.bioconductor.org/developers/how-to/unitTesting-guidelines/

context("bimodality")

test_that("bimodality scores work correctly", {

  expect_equal(bimodality_sarle(seq(100), type = "Sarle.finite.sample"), 0.2043795, tolerance=1.0e-3)
  
})

