#'@title Counts for all genes in bam data
#'
#'@description Function \code{countsForAllGenes} computes counts for all genes or gene parts in chosen gene data dictionary in the bam data set.
#'
#'@param data Bam data set 
#'@param geneData Data frame with positions of all genes and their parts, by default we use a \code{Araport11}.
#'@param genePart The part of gene we want to count apperances in data. By default we count apperances for whole gene.
#'@param range The distance behind the longest 3'UTR where we look for counts.
#'
#'@importFrom dplyr filter group_by %>%
#'
#'@export


countsForAllGenes <- function(data, geneData=directRNAExplorer::Araport11, genePart="gene", range=0){
  V3<-id<-gene<-id<-pos<-strand<-rname <- NULL
  geneInfo <- filter(geneData, V3==genePart)

  if(range>0){
    utrData <- filter(geneData, V3=="three_prime_UTR")
    utrData$length <- utrData$V5-utrData$V4
    utrData <- utrData%>%group_by(id)%>%filter(length==max(length)) 
    utrData <- utrData[!duplicated(utrData),]
  }
  
  countsgeneInfo <- data.frame(name = geneInfo$id, counts = 0)
  
  data <- data.frame(rname = data$rname,num=1:length(data$pos),pos=data$pos, strand=data$strand)
  
  
  for(i in 1:nrow(geneInfo)){
    chromosome <- geneInfo[i,"V1"]
    chromosome <- substr(chromosome, 4,5)
    
    data_chr <- subset(data, rname==chromosome)
    
    if(geneInfo[i,"V7"]=="+"){
      data_strand <- subset(data_chr, strand=="-")
    }else{
      data_strand <- subset(data_chr, strand=="+")
    }  
    if(range>0){
      utr_pos <- filter(geneInfo, id==geneInfo[i,"id"])
      utr_pos_e <- utr_pos$V5
      countsOutsideGene <- subset(data_strand, pos > utr_pos_e & pos <= (utr_pos_e + range))
      if(length(countsOutsideGene$pos)>=5){
      data_gen <-subset(data_strand, pos>=geneInfo[i, "V4"] & pos<=(utr_pos_e+range))
      }
    }
    if(range==0){
      data_gen <- subset(data_strand, pos>=geneInfo[i, "V4"] & pos<=geneInfo[i,"V5"])
    }
    
    counts_geneInfo[i,2] <- length(data_gen$pos)

  }
  
  return(counts_geneInfo)
}