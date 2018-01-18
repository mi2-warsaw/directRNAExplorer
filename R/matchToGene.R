#'@description Counts how many sequencing points were at known genes
#'
#'@title matchToGene
#'
#'@param positions Positions from sequencing.
#'@param start Gene start index.
#'@param stop Gene stop index.
#'@param geneName Gene name.
#'@param range How many nucleotide before \code{start} and after \code{stop} we include to genes.
#'
#'@examples
#'library(directRNAExplorer)
#'data <- dataChromosome1
#'dic <- TAIR10_genes
#'matchToGene(positions=unique(data$pos), start=dic$V4, stop=dic$V5, geneName=dic$id)
#'
#'@importFrom dplyr between
#'@export
matchToGene <- function(positions, start, stop, geneName, range=0){
  geneName <- geneName
  nrow <- length(geneName)
  counts <- as.data.frame(matrix(nrow=nrow,ncol=2))
  colnames(counts)<-c("count","name")
  if(range>0){
    start <- start-range
    stop <- stop+range
  }
  for (i in 1:length(start)){
    points <- positions[between(positions, start[i], stop[i])]
    counts$count[i]<-length(points)
  }
  counts$name <- geneName
  return(counts)
}
