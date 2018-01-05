context("Check input")

x <- c(1,2,3,4)

test_that("Invalid input",{
  expect_error(genesSummarize())
  expect_error(genesSummarize(x))
  
})

context("Check output")

data <- brmDataChromosome1[brmDataChromosome1$pos >2002000 & brmDataChromosome1$pos < 2018000,]

test_that("Output",{
  expect_is(genesSummarize(data), "data.frame")
  expect_length(genesSummarize(data),4)
})