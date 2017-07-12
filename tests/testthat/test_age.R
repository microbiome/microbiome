context("age operations")

test_that("Age grouping ok", {

  x <- c(1, 2, 10, 20, 50, NA, 100)	       
  
  expect_equal(length(unique(levels(group_age(x, "years")))), max(na.omit(x)) - min(na.omit(x)))
  expect_equal(length(unique(levels(group_age(x, "decades")))), 10)
  expect_equal(length(unique(levels(group_age(x, "even", 5)))), 5)    

})

