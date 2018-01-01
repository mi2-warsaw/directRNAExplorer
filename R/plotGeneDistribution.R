#'@title plotGeneDistribution
#'
#'@param gene character with gene name
#'@param bamDataFrame data frame converted to R using \code{bamToR()} function
#'@param geneData ata frame with positions of all genes and their names, by default we use a \code{TAIR10_genes}
#'@param range how many nucleotide before \code{start} and after \code{stop} we include to genes
#'
#'@importFrom dplyr filter
#'@importFrom tidyr separate
#'
#'@export

plotGeneDistribution <- function(gene, bamDataFrame, geneData=SequencingExplainer::TAIR10_genes, range=0){
  geneData <- filter(geneData, V3=="gene")
  geneData <- separate(geneData, col=V9, into=c("id","note","name"), sep="\\;")
  geneData$id <- substr(geneData$id, 4, length(geneData$id))
  geneInfo <- filter(geneData, id==gene)
  startPosition <- geneInfo$V4
  endPosition <- geneInfo$V5
  
  if(range>0){
    startPosition <- startPosition-range
    endPosition <- endPosition+range
  }
  
  chromosome <- geneInfo$V1
  chromosome <- factor(chromosome)
  chromosome <- substr(chromosome, 4,4)
  
  bam<- bamDataFrame@listData
  bam <- bam[c(1:11)]
  bam <- as.data.frame(bam)
  bamDataFrameGene <- filter(bam, pos>=startPosition & pos <=endPosition)
  bamDataFrameGene <- bamDataFrameGene[which(bamDataFrameGene$rname==chromosome),]
  
  plotCoverage(bamDataFrameGene, chromosome = chromosome)
  
}