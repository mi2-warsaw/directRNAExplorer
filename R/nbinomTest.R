#'Negative Binomial Test
#'
#'@param countsData Data frame containing counts for each 
#'@param condition Vector of conditions. In our data the type of mutation
#'@param ... Other arguments.
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom DESeq2 DESeq
#' @importFrom DESeq2 results
#'
#'@export


nbinomTest <- function(countsData, condition,...){
  data <- t(countsData)
  condition <- as.data.frame(condition)
  rownames(condition) <- colnames(data)
  
  datamx <- DESeqDataSetFromMatrix(countData = data,
                                   colData = condition,
                                   design = ~ condition)
  dds <- DESeq(datamx, ...)
  res <- results(dds)
  res2 <- as.data.frame(res@listData)
  rownames(res2) <- res@rownames
  res2$id <- rownames(res2)
  return(res2)
}