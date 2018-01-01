#'@description Function genesSumarize counts the coverage of each gene
#'
#'@title genesSumarize
#'
#'@param bamDataFrame data frame converted to R using \code{bamToR()} function
#'@param geneData data frame with positions of all genes and their names, by default we use a \code{TAIR10_genes}
#'@param chromosome optional, number of chromosome
#'@param range how many nucleotide before \code{start} and after \code{stop} we include to genes
#'
#'@importFrom dplyr filter
#'@importFrom tidyr separate
#'
#'@export

genesSummarize <- function(bamDataFrame, geneData=SequencingExplainer::TAIR10_genes, chromosome=NULL, range=0){
  if(!is.null(chromosome)){
    geneData <- filter(geneData, V1==paste("Chr",chromosome))
  }
  geneData <- filter(geneData, V3=="gene")
  geneData <- separate(geneData, col=V9, into=c("id","note","name"), sep="\\;")
  geneData$id <- substr(geneData$id, 4, length(geneData$id))
  counts <- matchToGene(positions = bamDataFrame$pos, start=geneData$V4, stop=geneData$V5, geneName=geneData$id, range = range)
  uniqueCounts <- matchToGene(positions = unique(bamDataFrame$pos),start=geneData$V4, stop=geneData$V5, geneName=geneData$id, range = range)
  genesSummary <- data.frame(name = geneData$id, counts = counts$count, uniqueCounts = uniqueCounts$count, length = (geneData$V5 - geneData$V4))

  return(genesSummary)
}
