context("taxonomic map")

test_that("levelmap works correctly", {

  data(dietswap)
  tt <- tax_table(dietswap)

  # For phyloseq

  expect_equal(map_levels('Akkermansia', from = 'Genus', to = 'Phylum', tt), "Verrucomicrobia")
  expect_equal(map_levels('Verrucomicrobia', from = 'Phylum', to = 'Genus', tt)$Verrucomicrobia, "Akkermansia")
})

