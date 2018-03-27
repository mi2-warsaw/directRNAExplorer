#'@description Function genesSumarize counts the coverage of each gene
#'
#'@title genesSumarize
#'
#'@param bamDataFrame Data frame converted to R using \code{bamToR()} function
#'@param geneData Data frame with positions of all genes and their names, by default we use a \code{TAIR10_genes_tidy}
#'@param chromosome Optional, number of chromosome
#'@param range How many nucleotide before \code{start} and after \code{stop} we include to genes.
#'@param strand On which strand we want to compute statistics, by default we use the \code{sense} strand, we can also compute statistics on the \code{antisense} strand.
#'@param genePart Part of gene on which we compute the summary
#'
#'@examples
#'library(directRNAExplorer)
#'data <- brmDataChromosome1[brmDataChromosome1$pos >2002000 & brmDataChromosome1$pos < 2018000,]
#'genesSummarize(data)
#'
#'@importFrom dplyr filter
#'@importFrom tidyr separate
#'
#'@export

genesSummarize <- function(bamDataFrame, geneData=directRNAExplorer::TAIR10_genes_tidy, chromosome=NULL, range=0, strand="sense", genePart="gene"){
  V1 <- NULL
  if(!is.null(chromosome)){
    geneData <- filter(geneData, V1==paste("Chr",chromosome, sep=""))
  }
  
  counts <- matchToGene(positions = data.frame(pos = bamDataFrame$pos, strand = bamDataFrame$strand), start=geneData$V4, stop=geneData$V5, geneName=geneData$id, strandData=geneData$V7, strand = strand, range = range)
  uniqueCounts <- matchToGene(positions = unique(data.frame(pos = bamDataFrame$pos, strand = bamDataFrame$strand)),start=geneData$V4, stop=geneData$V5, geneName=geneData$id, strandData=geneData$V7, strand = strand,range = range)
  
  genesSummary <- data.frame(name = geneData$id, counts = counts$count, uniqueCounts = uniqueCounts$count, length = (geneData$V5 - geneData$V4))

  return(genesSummary)
}
