#'@title T-student test
#'
#'@description Function \code{tTest} computes statistics for genes with division into two groups of probes. 
#'
#'@param data Data set containing counts for each gene, columns corresponds to genes, rows to samples
#'@param condition A two - levels vector of condition values 
#'@param ... Other arguments
#'
#'@importFrom limma lmFit  eBayes contrasts.fit topTable makeContrasts
#'
#'@export

tTest <- function(data, condition, ...) {
  data <- t(data)
  levels <- unique(condition)
  design <- data.frame(c1 = ifelse(condition == levels[1], 1, 0),
                       c2 = ifelse(condition == levels[2], 1, 0))
  colnames(design) <- levels
  fit <- lmFit(data, design)
  forms <- paste0(colnames(design)[1], "-", colnames(design)[2])
  contMatrix <- makeContrasts(contrasts = forms, levels = design)
  fit2 <- contrasts.fit(fit, contMatrix)
  fit2 <- eBayes(fit2)
  table <- topTable(fit2, number = Inf, coef = 1, ...)
  table$id <- rownames(table)
  colnames(table) <- c("log2fold", "mean", "t.stat", "pval", "padj", "B", "id")
  table <- table[, c(7, 2, 1, 4, 5)]
  return(table)
}