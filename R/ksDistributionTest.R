#'@title ksDistributionTest
#'
#'@description Kolmogorov-Smirnov test for gene from two datasets
#'
#'@param data1 First data set which contains informations about selected gene.
#'@param data2 Second data set which contains informations about selected gene.
#'@param gene Gene name.
#'@param geneData Data frame with positions of all genes and their names, by default we use a \code{TAIR10_genes}
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


ksDistributionTest <- function(data1, data2, gene, geneData = directRNAExplorer::TAIR10_genes, genePart="gene"){
  id<-V4<-V5<-V3 <-V7<- NULL
  geneInfo <- filter(geneData, id==gene)
  geneInfo <- filter(geneInfo, V3==genePart)
  if(dim(geneInfo)[1]>1){
    startPosition <-min(geneInfo$V4)
    endPosition <- max(geneInfo$V5)
  }else{
    startPosition <- geneInfo$V4
    endPosition <- geneInfo$V5
  }
  
  strand <- geneInfo$V7
  geneRange <- c(startPosition, endPosition)
 # if((geneRange[1]<min(data1$pos) || geneRange[2]>max(data1$pos)) || (geneRange[1]<min(data2$pos) || geneRange[2]>max(data2$pos))){
#    stop("Your data sets don't contain selected gene")
 # }
  
  data1 <- data1[data1$pos>=geneRange[1] & data1$pos<=geneRange[2],]
  data2 <- data2[data2$pos>=geneRange[1] & data2$pos<=geneRange[2],]
  
  if(strand=="+"){
    position <- c(data1$pos, data2$pos)
    strand <- c(as.character(data1$strand), as.character(data2$strand))
    condition <- c(rep("type1", length(data1$pos)), rep("type2", length(data2$pos)))
    
    dt <- data.frame(position, strand, condition)
    dt_pos <- dplyr::filter(dt, strand == "-")
    if(length(dt_pos[dt_pos$condition=="type1","position"])<2 || length(dt_pos[dt_pos$condition=="type2", "position"])<2){
      return(data.frame(p.value=NA,statistic=NA))
    }else{
  return(ks.test(dt_pos[dt_pos$condition=="type1","position"], dt_pos[dt_pos$condition=="type2", "position"]))  
    }
  }
  if(strand=="-"){
    position <- c(data1$pos, data2$pos)
    strand <- c(as.character(data1$strand), as.character(data2$strand))
    condition <- c(rep("type1", length(data1$pos)), rep("type2", length(data2$pos)))
    
    dt <- data.frame(position, strand, condition)
    dt_neg <- dplyr::filter(dt, strand == "+")
    if(length(dt_neg[dt_neg$condition=="type1","position"])<2 || length(dt_neg[dt_neg$condition=="type2", "position"])<2){
      return(data.frame(p.value=NA,statistic=NA))
    }else{
    return(ks.test(dt_neg[dt_neg$condition=="type1","position"], dt_neg[dt_neg$condition=="type2", "position"]))
  }
  } 
}