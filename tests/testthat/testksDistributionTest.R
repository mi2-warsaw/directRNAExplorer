databrm <- brmDataChromosome1[brmDataChromosome1$pos >2002000 & brmDataChromosome1$pos < 2018000,]
datawt <- dataChromosome1[dataChromosome1$pos >2002000 & dataChromosome1$pos < 2018000,]

context("Check input")

test_that("Invalid input",{
  expect_error(ksDistributionTest())
  expect_error(ksDistributionTest(databrm, datawt))
  expect_error(ksDistributionTest(databrm, gene = "AT1G06560"))
  
})

context("Check output")

test_that("Output",{
  expect_is(ksDistributionTest(databrm, datawt, gene = "AT1G06560"), "htest")
  expect_length(ksDistributionTest(databrm, datawt, gene = "AT1G06560"),5)
  expect_is(ksDistributionTest(databrm, datawt, gene = "AT1G06560", genePart = "three_prime_UTR", strand = "neg"),"htest")
  expect_is(ksDistributionTest(databrm, datawt, gene = "AT1G06560", genePart = "three_prime_UTR"),"htest")
  expect_length(ksDistributionTest(databrm, datawt, gene = "AT1G06560", genePart = "three_prime_UTR"),5)
})