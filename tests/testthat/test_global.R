context("alpha")

test_that("alpha indices work correctly", {

  data(peerj32)
  pseq <- peerj32$phyloseq
  
  # expect_true(is.data.frame(alpha(pseq)))
  expect_true(is.data.frame(rarity(pseq)))
  expect_true(is.data.frame(dominance(pseq)))
  expect_true(is.data.frame(alpha(pseq)))
  expect_true(is.data.frame(evenness(pseq)))

  expect_true(all(rownames(abundances(pseq)) == taxa(pseq)))

  expect_true(all(unlist(richness(pseq, index='observed') ) == unlist(alpha(pseq, index='observed') )))      

  expect_true(max(low_abundance(transform(pseq, "compositional"))) <= 1)

  expect_true(all(c("evenness_simpson", "dominance_simpson") %in%
      names(alpha(pseq, index = "simpson"))))
  expect_error(suppressWarnings(alpha(transform(pseq, "clr"), index = "fisher")))

})





test_that("dominance works correctly", {

  data(peerj32)
  pseq <- peerj32$phyloseq
  
  expect_true(is.null(dominance(pseq, index = "shannon")))
  
  expect_true(all(dominance(pseq)$dbp == dominance(pseq, index = "DBP")))  

  expect_true(all(dominance(pseq, index = NULL, relative = FALSE, rank = 1) ==
               dominance(pseq, index = "absolute")))

  expect_true(all(dominance(pseq, index = NULL, relative = TRUE, rank = 1) ==
               dominance(pseq, index = "relative")))

  expect_true(all(dominance(pseq, index = NULL, relative = TRUE, rank = 1) == 
               dominance(pseq, index = "DBP")))

  expect_true(all(dominance(pseq, index = NULL, relative = TRUE, rank = 2) ==
               dominance(pseq, index = "DMN")))


})



test_that("evenness works correctly", {

  data(peerj32)
  pseq <- peerj32$phyloseq  
  expect_true(all(evenness(abundances(pseq)) == evenness(pseq)))

})



test_that("alpha indices are harmonized", {

  data(peerj32)
  pseq <- peerj32$phyloseq
  expect_true(is.data.frame(evenness(abundances(pseq), index = "camargo")))
  expect_true(is.data.frame(alpha(abundances(pseq), index = "shannon")))

})





test_that("richness works correctly", {

  data(peerj32)
  pseq <- peerj32$phyloseq  
  expect_true(all(richness(abundances(pseq)) == richness(pseq)))
  expect_true(all(richness(pseq, detection = 0)[,1] == richness(pseq)[,1])  )


})


