countsPerGene <- function(datawt = list(), databrm = list(), datac = list(),  gene, geneData=directRNAExplorer::TAIR10_genes, genePart = "gene", type="unique"){
  pos<-id<-V4<-V5<-V3 <-name<- NULL
  geneInfo <- filter(geneData, id==gene)
  geneInfo <- filter(geneInfo, V3==genePart)
  if(dim(geneInfo)[1]>1){
    startPosition <-min(geneInfo$V4)
    endPosition <- max(geneInfo$V5)
  }else{
    startPosition <- geneInfo$V4
    endPosition <- geneInfo$V5
  }
  
  
  chromosome <- geneInfo$V1
  chromosome <- factor(chromosome)
  chromosome <- substr(chromosome, 4,4)
  chromosome <- chromosome[1]
  ####
  countsData1 <- data.frame(counts=0, condition = rep("typewt", length(datawt)), name = names(datawt))
  
  
  countsData2 <- data.frame(counts=0, condition = rep("typebrm", length(databrm)), name = names(databrm))
 
  
  countsData3 <- data.frame(counts=0, condition = rep("typec", length(datac)), name = names(datac))
  
  
  
  for(i in 1:length(datawt)){
    datawt[[i]] <- subset(datawt[[i]], pos>=startPosition & pos<=endPosition)
  }
  
  for(i in 1:length(databrm)){
    databrm[[i]] <- subset(databrm[[i]], pos>=startPosition & pos<=endPosition)
  }
  
  for(i in 1:length(datac)){
    datac[[i]] <- subset(datac[[i]], pos>=startPosition & pos<=endPosition)
  }
  
  
  for(i in 1:length(datawt)){
    summaryData <- genesSummarize(datawt[[i]])
    summaryData <- subset(summaryData, name==gene)
    if(type=="unique"){
      countsData1[i,1] <- summaryData[1,3]
    }else{
      countsData1[i,1] <- summaryData[1,2]
    }
  }
  
  for(i in 1:length(databrm)){
    summaryData <- genesSummarize(databrm[[i]])
    summaryData <- subset(summaryData, name==gene)
    if(type=="unique"){
      countsData2[i,1] <- summaryData[1,3]
    }else{
      countsData2[i,1] <- summaryData[1,2]
    }
  }
  
  for(i in 1:length(datac)){
    summaryData <- genesSummarize(datac[[i]])
    summaryData <- subset(summaryData, name==gene)
    if(type=="unique"){
      countsData3[i,1] <- summaryData[1,3]
    }else{
      countsData3[i,1] <- summaryData[1,2]
    }
  }
  
  return(rbind(countsData1, countsData2, countsData3))
}