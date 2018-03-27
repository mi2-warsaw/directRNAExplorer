#'@description Counts how many sequencing points were at known genes
#'
#'@title matchToGene
#'
#'@param positions Positions from sequencing and informations about strand.
#'@param start Gene start index.
#'@param stop Gene stop index.
#'@param geneName Gene name.
#'@param range How many nucleotide before \code{start} and after \code{stop} we include to genes.
#'@param strandData Vector with information that gene is on \code{sense} or \code{antisense} strand.
#'@param strand On which strand we want to compute statistics, we can compute statistics on the \code{sense} or \code{antisense} strand.
#'
#'@examples
#'\dontrun{
#'library(directRNAExplorer)
#'data <- dataChromosome1
#'dic <- TAIR10_genes
#'matchToGene(positions=unique(data$pos), start=dic$V4, stop=dic$V5, geneName=dic$id)
#'}
#'@importFrom dplyr between filter
#'@export
matchToGene <- function(positions, start, stop, geneName, strandData, strand="sense", range=0){
  geneName <- geneName
  nrow <- length(geneName)
  counts <- as.data.frame(matrix(nrow=nrow,ncol=2))
  colnames(counts)<-c("count","name")
  if(range>0){
    start <- start-range
    stop <- stop+range
  }
  for (i in 1:length(start)){
    if(strand=="sense"){
      if(strandData[i]=="+"){
        positions <- subset(positions, strand=="-")
      }else{
        positions <- subset(positions, strand=="+")
      }
      
    }
    if(strand=="antisense"){
      if(strandData[i]=="+"){
        positions <- subset(positions, strand=="+")
      }else{
        positions <- subset(positions, strand=="-")
      }
    }
    points <- positions[,1][between(positions[,1], start[i], stop[i])]
    counts$count[i]<-length(points)
  }
  counts$name <- geneName
  return(counts)
}
