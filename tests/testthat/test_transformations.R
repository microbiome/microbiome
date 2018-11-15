context("phyloseq transformations")

test_that("transform works correctly", {
  # Hellinger distance is sqrt of the relative abundance
  # Note the different scales (100x) of these two measures.
  
  data(dietswap)
  
  hel <- otu_table(transform(dietswap, "hellinger"))^2
  rel <- otu_table(transform(dietswap, "compositional"))
  expect_equal(max(abs(hel - rel)), 0, tolerance = 1e-15)

  expect_equal(ntaxa(transform(dietswap, "identity")), ntaxa(dietswap))

  z <- transform(transform(dietswap, "shift", shift = 1), "Z")
  expect_equal(ntaxa(z), ntaxa(dietswap))
  expect_equal(sum(apply(abundances(z), 1, mean)), 0, tolerance = 1e-10)
  
  expect_equal(ntaxa(transform(dietswap, "clr")), ntaxa(dietswap))
  expect_equal(ntaxa(transform(transform(dietswap, "shift", shift = 1), "log10")), ntaxa(dietswap))
  expect_true(all(transform(abundances(transform(dietswap, "shift", shift = 1))[1:3, 1:3], "log10") == transform(abundances(dietswap)[1:3, 1:3], "log10p")))    

  expect_equal(ntaxa(transform(dietswap, "shift", shift = 1)), ntaxa(dietswap))      
  expect_equal(ntaxa(transform(dietswap, "compositional")), ntaxa(dietswap))
  expect_true(sum(colSums(abundances(transform(dietswap, "compositional"))) - 1) < 1e-15)

})

