context("global")

test_that("global indices work correctly", {

  data(peerj32)
  pseq <- peerj32$phyloseq
  
  expect_true(is.vector(global(pseq, measures = "Observed", split = FALSE)))
  expect_true(is.data.frame(global(pseq, split = TRUE)))
  
  expect_true(is.vector(rarity(pseq)))
 
  expect_true(max(low_abundance(transform(pseq, "compositional"))) <= 1)  

  expect_equal(rare_abundance(pseq, detection = 0.1/100, prevalence = 50/100),
               1 - core_abundance(pseq, detection = 0.1/100, prevalence = 50/100))

})




test_that("dominance works correctly", {

  data(peerj32)
  pseq <- peerj32$phyloseq
  
  expect_equal(dominance(pseq),
               dominance(pseq, index = "DBP"))

  expect_equal(dominance(pseq, index = NULL, relative = FALSE, rank = 1),
               dominance(pseq, index = "absolute_dominance"))

  expect_equal(dominance(pseq, index = NULL, relative = TRUE, rank = 1),
               dominance(pseq, index = "relative_dominance"))

  expect_equal(dominance(pseq, index = NULL, relative = TRUE, rank = 1),
               dominance(pseq, index = "DBP"))

  expect_equal(dominance(pseq, index = NULL, relative = TRUE, rank = 2, aggregate = TRUE),
               dominance(pseq, index = "DMN"))


})


