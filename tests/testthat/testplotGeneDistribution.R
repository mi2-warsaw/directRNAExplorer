
data <- brmDataChromosome1[brmDataChromosome1$pos >2002000 & brmDataChromosome1$pos < 2018000,]

context("Check input")

test_that("Invalid input",{
  expect_error(plotGeneDistribution())
  expect_error(plotGeneDistribution(data))
  
})

context("Check output")

test_that("Output",{
  expect_is(plotGeneDistribution("AT1G06560", data), "gtable")
  expect_is(plotGeneDistribution("AT1G06560", data, genePart="three_prime_UTR"), "gtable")
  expect_is(plotGeneDistribution("AT1G06560", data, genePart="three_prime_UTR", range = 10, type = "density"), "gtable")
  expect_length(plotGeneDistribution("AT1G06560", data), 34)
})