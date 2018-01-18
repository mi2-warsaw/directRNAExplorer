countsInGene <- function(data, geneStart, geneStop, strand){
  if(strand == "pos"){
    countsData <- data[apply(as.data.frame(data$flag),1, check_pos), 'pos']
  }
  if(strand=="neg"){
    countsData <- data[apply(as.data.frame(data$flag),1, check_neg), 'pos']
  }
  
  countsData <- as.data.frame(table(countsData))
  colnames(countsData)[1] <- "pos"
  
  #if(dim(countsData)[1]==0){
  #  stop(paste("Selected gene has no counts on the", strand, "strand."))
  #}
  
  if(dim(countsData)[1]==0){
    countsData <- data.frame(pos = factor(c(geneStart:geneStop)), Freq = 0)
    
  }
  countsStrandData <- data.frame(pos=factor(c(geneStart:geneStop)))
  countsStrandData <- dplyr::left_join(countsStrandData, countsData, by="pos")
  countsStrandData[is.na(countsStrandData)]<-0
  return(countsStrandData)
}