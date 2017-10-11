context("BMI operations")

test_that("BMI grouping ok", {

  x <- 15:50
  
  expect_equal(length(unique(levels(group_bmi(x, "standard")))), 7)
  expect_true(group_bmi(x, "standard")[[10]] == "lean")  
  expect_equal(length(unique(levels(group_bmi(x, "even", 5)))), 5)

  expect_true(group_bmi(18, "standard") == "underweight")
  expect_true(group_bmi(18.5, "standard") == "lean")    
  expect_true(group_bmi(20, "standard") == "lean")
  expect_true(group_bmi(22, "standard") == "lean")  
  expect_true(group_bmi(25, "standard") == "overweight")
  expect_true(group_bmi(28, "standard") == "overweight")  
  expect_true(group_bmi(30, "standard") == "obese")
  expect_true(group_bmi(32, "standard") == "obese")  
  expect_true(group_bmi(35, "standard") == "severe")
  expect_true(group_bmi(37, "standard") == "severe")
  expect_true(group_bmi(40, "standard") == "morbid")
  expect_true(group_bmi(42, "standard") == "morbid")  
  expect_true(group_bmi(45, "standard") == "super")
  expect_true(group_bmi(47, "standard") == "super")  
  expect_true(group_bmi(50, "standard") == "super")
  expect_true(group_bmi(52, "standard") == "super")  
  expect_true(group_bmi(55, "standard") == "super")      

})

