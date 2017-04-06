context("diversity")

test_that("diversity works correctly", {

  data(peerj32)
  pseq <- peerj32$phyloseq
  
  expect_true(is.vector(diversity(pseq, split = FALSE)))
  expect_true(is.data.frame(diversity(pseq, split = TRUE)))
  
  expect_true(is.vector(rarity(pseq, detection = 0.1/100, prevalence = 50/100, split = TRUE)))
  expect_true(length(rarity(pseq, detection = 0.1/100, prevalence = 50/100, split = FALSE)) == 0)
  
  expect_true(min(dominance(pseq)) >= 0)
  expect_true(max(dominance(pseq)) <= 1)

  expect_true(is.integer(dominance(pseq, threshold = 0.5)[[1]]))

  expect_equal(rarity(pseq, detection = 0.1/100, prevalence = 50/100),
               1 - core_abundance(pseq, detection = 0.1/100, prevalence = 50/100))  

})


