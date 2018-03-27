#'pvalTable
#'
#'Create data.frame with pvalues and value of statistic for chosen statistical test
#'
#'@param data1 First data set which contains informations about selected gene.
#'@param data2 Second data set which contains informations about selected gene.
#'@param chromosome Chromosome name.
#'@param geneData Data frame with positions of all genes and their names, by default we use a \code{TAIR10_genes}.
#'@param type Type of chosen test
#'@param genePart The part of gene we want to visualize.
#'@param ... Optional arguments
#'
#'@export


pvalTable <- function(data1, data2, chromosome, geneData=directRNAExplorer::TAIR10_genes,  type="KS", genePart = "gene", ...){
  V1<-V3<-V4<-V5<-rname <- NULL
  geneData$V1 <- factor(substr(geneData$V1, 4,4))
  geneData <- dplyr::filter(geneData, V1==chromosome)
  if(genePart!="gene"){
    geneData <- dplyr::filter(geneData, V3 == genePart) 
  }else{
  geneData <- dplyr::filter(geneData, V3 == "gene")
  }
  rangeData <- c(max(min(data1$pos), min(data2$pos)), min(max(data1$pos), max(data2$pos)))
  geneData2 <- dplyr::filter(geneData, V4>rangeData[1] & V5<rangeData[2])
  
  genes <- unique(geneData2$id)
  data1 <- data1[data1$rname==chromosome, ]
  data2 <- data2[data2$rname== chromosome, ]
  
  result <- data.frame(gene = genes, pval=0, statistic=0)
  if(type=="KS"){
    for (i in 1:length(genes)){
      test <- ksDistributionTest(data1,data2,result[i,1], genePart = genePart)
      result[i,2] <- test$p.value
      result[i,3] <- test$statistic
    }
  }
  
  return(result)
  
}