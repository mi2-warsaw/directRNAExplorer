context("check input")

x <- c(1,2,3,4)


test_that("Invalid argument",{
  expect_error(bamToR())
  expect_error(bamToR(x))

})
