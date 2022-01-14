context("alpha")

test_that("alpha indices work correctly", {

  data(dietswap)
  pseq <- dietswap

  expect_equal(ntaxa(aggregate_rare(dietswap, level = 'Genus', detection = 0.1/100, prevalence = 5/100)), 116)

  # expect_true()

  expect_true(is.vector(core_abundance(pseq, detection = 0.1/100, prevalence = 50/100)))
  
  expect_equal(length(core_abundance(prune_samples("Sample-1", pseq), detection = 10/100, prevalence = 50/100)), 1)


})
