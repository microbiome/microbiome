context("phyloseq transformations")

test_that("transform_phyloseq works correctly", {
  # Hellinger distance is sqrt of the relative abundance
  # Note the different scales (100x) of these two measures.
  data(atlas1006)
  hel <- otu_table(transform_phyloseq(atlas1006, "hellinger"))^2
  rel <- otu_table(transform_phyloseq(atlas1006, "relative.abundance"))/100
  expect_equal(max(abs(hel - rel)), 0, tolerance = 1e-15)
})

