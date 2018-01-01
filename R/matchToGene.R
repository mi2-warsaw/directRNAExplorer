#'@description Counts how many sequencing points were at known genes
#'
#'@title matchToGene
#'
#'@param positions positions from sequencing
#'@param start gene start index
#'@param stop gene stop index
#'@param geneName gene name
#'@param range how many nucleotide before \code{start} and after \code{stop} we include to genes
#'
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
