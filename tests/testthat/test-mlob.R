library(testthat)

# Test that mlob runs correctly and returns proper object
test_that("mlob runs without error and returns correct class", {
  result <- mlob(Sepal.Length ~ Sepal.Width + Petal.Length, 
                 data = iris, group = "Species")
  expect_s3_class(result, "mlob_result")
  expect_true(is.list(result))
})

# Test expected components of the result
test_that("mlob_result contains expected components", {
  result <- mlob(Sepal.Length ~ Sepal.Width + Petal.Length, 
                 data = iris, group = "Species")
  expect_true("beta_b" %in% names(result$Coefficients))
})

# Test summary method
test_that("summary.mlob_result produces output", {
  result <- mlob(Sepal.Length ~ Sepal.Width + Petal.Length, 
                 data = iris, group = "Species")
  expect_output(summary(result), "Summary of Coefficients")
})

# Test print method
test_that("print.mlob_result produces output", {
  result <- mlob(weight ~ Time, data = ChickWeight, group = "Chick")
  expect_output(print(result), "Call:")
})


# Test with jackknife = TRUE
test_that("mlob works with jackknife = TRUE", {
  result <- mlob(Sepal.Length ~ Sepal.Width + Petal.Length, 
                 data = iris, group = "Species", jackknife = TRUE)
  expect_s3_class(result, "mlob_result")
})

# Test behavior with custom balancing.limit and punish.coeff
test_that("mlob accepts custom balancing.limit and punish.coeff", {
  result <- mlob(Sepal.Length ~ Sepal.Width + Petal.Length, 
                 data = iris, group = "Species", balancing.limit = 0.3, punish.coeff = 1.5)
  expect_s3_class(result, "mlob_result")
})

# Test handling of invalid group variable
test_that("mlob throws error on invalid group variable", {
  expect_error(
    mlob(Sepal.Length ~ Sepal.Width + Petal.Length, 
         data = iris, group = "NotARealGroup")
  )
})

# Test formula input validation
test_that("mlob throws error when formula is invalid", {
  expect_error(
    mlob("not a formula", data = iris, group = "Species")
  )
})
