#'@title ksDistributionTest
#'
#'@description Kolmogorov-Smirnov test for gene from two datasets
#'
#'@param data1 First data set which contains informations about selected gene.
#'@param data2 Second data set which contains informations about selected gene.
#'@param gene Gene name
#'@param geneData Data frame with positions of all genes and their names, by default we use a \code{TAIR10_genes_tidy}
#'@param strand On which strand we want to compute test
#'
#'
#'@importFrom dplyr filter
#'@importFrom stats ks.test
#'
#'@export


ksDistributionTest <- function(data1, data2, gene, geneData = SequencingExplainer::TAIR10_genes_tidy, strand="both"){
  id<-V4<-V5 <- NULL
  geneData <- filter(geneData, id == gene)
  geneRange <- c(geneData$V4, geneData$V5)
  if((geneRange[1]<min(data1$pos) || geneRange[2]>max(data1$pos)) || (geneRange[1]<min(data2$pos) || geneRange[2]>max(data2$pos))){
    stop("Your data sets don't contain selected gene")
  }
  
  #data1 <- subset(data1, pos>=geneRange[1] & pos<=geneRange[2])
  #data2 <- subset(data2, pos>=geneRange[1] & pos<=geneRange[2])
  
  data1 <- data1[data1$pos>=geneRange[1] & data1$pos<=geneRange[2],]
  data2 <- data2[data2$pos>=geneRange[1] & data2$pos<=geneRange[2],]
  
  if(strand=="pos"){
    posStrandData1 <- countsInGene(data1, geneRange[1], geneRange[2], strand="pos")
    posStrandData2 <- countsInGene(data2, geneRange[1], geneRange[2], strand="pos")
  return(ks.test(posStrandData1$Freq, posStrandData2$Freq))
  }
  
  if(strand=="neg"){
    negStrandData1 <- countsInGene(data1, geneRange[1], geneRange[2], strand="neg")
    negStrandData2 <- countsInGene(data2, geneRange[1], geneRange[2], strand="neg")
    return(ks.test(negStrandData1$Freq, negStrandData2$Freq))
  }
  
  if(strand=="both"){
    posStrandData1 <- countsInGene(data1, geneRange[1], geneRange[2], strand="pos")
    posStrandData2 <- countsInGene(data2, geneRange[1], geneRange[2], strand="pos")
    negStrandData1 <- countsInGene(data1, geneRange[1], geneRange[2], strand="neg")
    negStrandData2 <- countsInGene(data2, geneRange[1], geneRange[2], strand="neg")
    return(list(ks.test(posStrandData1$Freq, posStrandData2$Freq), ks.test(negStrandData1$Freq, negStrandData2$Freq)))
  }else{
    stop("You must set the strand argument")
  }
}