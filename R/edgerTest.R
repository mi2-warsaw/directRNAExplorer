#'EdgeR based test for differential expression analysis
#'
#'@description Function \code{test_edger} computes the statistics for quasi-likelihood F-tests or likelihood ratio tests. 
#'
#'@param countsData Data frame containing counts for each 
#'@param condition Vector of conditions. In our data the type of mutation
#'@param type Type of test "lrt" for likelihood ratio tests, "qlf" for quasi-likelihood F-tests.
#'@param ... Other arguments.
#'
#'@importFrom edgeR DGEList calcNormFactors estimateDisp glmFit glmLRT glmQLFit glmQLFTest
#'@importFrom stats model.matrix
#'
#'@export

edgerTest <- function(countsData, condition, type="qlf",...){
  y <- DGEList(counts = t(countsData), group = condition)
  y <- calcNormFactors(y)
  design <- model.matrix(~condition)
  y <- estimateDisp(y, design)
  if (type == "lrt") {
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, coef = 2)
    result <- lrt@.Data[[14]]
  }
  if (type == "qlf") {
    fit <- glmQLFit(y, design)
    qlf <- glmQLFTest(fit, coef = 2)
    result <- qlf@.Data[[17]]
  }
  result$id <- rownames(result)
    colnames(result)[4] <- "pvalue"
    colnames(result)[1] <- "log2FoldChange"
    return(result)
}