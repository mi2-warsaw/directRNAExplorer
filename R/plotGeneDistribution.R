#'@title plotGeneDistribution
#'
#'@description Plot of distribution on selected gene
#'
#'@param gene Character with gene name.
#'@param bamDataFrame Data frame converted to R using \code{bamToR()} function.
#'@param geneData Data frame with positions of all genes and their names, by default we use a \code{TAIR10_genes}.
#'@param range How many nucleotide before \code{start} and after \code{stop} we include to genes.
#'@param genePart The part of gene we want to visualize.
#'@param ... Optional arguments.
#'
#'@importFrom dplyr filter
#'@importFrom tidyr separate
#'
#'@examples
#'library(SequencingExplainer)
#'data <- brmDataChromosome1[brmDataChromosome1$pos >2002000 & brmDataChromosome1$pos < 2018000,]
#'plotGeneDistribution("AT1G06560", data, genePart="three_prime_UTR")
#'
#'@export

plotGeneDistribution <- function(gene, bamDataFrame, geneData=SequencingExplainer::TAIR10_genes, range=0, genePart="gene", ...){
  pos<-id<-V4<-V5<-V3 <- NULL
  geneInfo <- filter(geneData, id==gene)
  geneInfo <- filter(geneInfo, V3==genePart)
  if(dim(geneInfo)[1]>1){
    startPosition <-min(geneInfo$V4)
    endPosition <- max(geneInfo$V5)
  }else{
    startPosition <- geneInfo$V4
    endPosition <- geneInfo$V5
  }
  
  if(range>0){
    startPosition <- startPosition-range
    endPosition <- endPosition+range
  }
  
  chromosome <- geneInfo$V1
  chromosome <- factor(chromosome)
  chromosome <- substr(chromosome, 4,4)
  chromosome <- chromosome[1]
  
  bam<- bamDataFrame@listData
  bam <- bam[c(1:11)]
  bam <- as.data.frame(bam)
  bamDataFrameGene <- filter(bam, pos>=startPosition & pos <=endPosition)
  bamDataFrameGene <- bamDataFrameGene[which(bamDataFrameGene$rname==chromosome),]
  
  return(plotCoverage(bamDataFrameGene, chromosome = chromosome, start=startPosition, stop=endPosition))
  
}