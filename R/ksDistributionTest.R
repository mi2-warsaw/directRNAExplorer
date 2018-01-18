#'@title ksDistributionTest
#'
#'@description Kolmogorov-Smirnov test for gene from two datasets
#'
#'@param data1 First data set which contains informations about selected gene.
#'@param data2 Second data set which contains informations about selected gene.
#'@param gene Gene name.
#'@param geneData Data frame with positions of all genes and their names, by default we use a \code{TAIR10_genes}.
#'@param strand On which strand we want to compute test.
#'@param genePart The part of gene we want to test. By default we conduct test for whole gene.
#'
#'@examples
#'library(directRNAExplorer)
#'databrm <- brmDataChromosome1[brmDataChromosome1$pos >2002000 & brmDataChromosome1$pos < 2018000,]
#'datawt <- dataChromosome1[dataChromosome1$pos >2002000 & dataChromosome1$pos < 2018000,]
#'ksDistributionTest(databrm, datawt, gene = "AT1G06560")
#'
#'@importFrom dplyr filter
#'@importFrom stats ks.test
#'
#'@export


ksDistributionTest <- function(data1, data2, gene, geneData = SequencingExplainer::TAIR10_genes, strand="pos", genePart="gene"){
  id<-V4<-V5<-V3 <- NULL
  geneInfo <- filter(geneData, id==gene)
  geneInfo <- filter(geneInfo, V3==genePart)
  if(dim(geneInfo)[1]>1){
    startPosition <-min(geneInfo$V4)
    endPosition <- max(geneInfo$V5)
  }else{
    startPosition <- geneInfo$V4
    endPosition <- geneInfo$V5
  }
  
  geneRange <- c(startPosition, endPosition)
 # if((geneRange[1]<min(data1$pos) || geneRange[2]>max(data1$pos)) || (geneRange[1]<min(data2$pos) || geneRange[2]>max(data2$pos))){
#    stop("Your data sets don't contain selected gene")
 # }
  
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