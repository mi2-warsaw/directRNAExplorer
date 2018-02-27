#'@description Histogram plot for chosen part of RNA sequence in two datasets
#'
#'@title plotHistogram
#'
#'@param gene Character with gene name.
#'@param data1 Data frame converted to R using \code{bamToR()} function.
#'@param geneData Data frame with positions of all genes and their names, by default we use a \code{TAIR10_genes}.
#'@param range How many nucleotide before \code{start} and after \code{stop} we include to genes.
#'@param genePart The part of gene we want to visualize.
#'
#'@importFrom ggplot2 ggplot theme_bw aes theme element_blank labs ylab scale_y_reverse xlim geom_histogram
#'@export

plotHistogram <- function(gene, data1, geneData=directRNAExplorer::TAIR10_genes, range=0, genePart="gene"){
  pos<-id<-V4<-V5<-V3<-position<-type <- NULL
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
  
  if(geneInfo[,"V7"]=="+"){
    data1 <- data1[which(data1$strand=="-"),]
  }else{
    data1 <- data1[which(data1$strand=="+"),]
  }
  
  data1 <- data1[which(data1$pos>=startPosition & data1$pos <=endPosition),]
  
  position <- data1$pos
  strand <- as.character(data1$strand)
  
  data <- data.frame(position, strand)
  
  range <- c(startPosition, endPosition)
  
  
  plot1 <- ggplot(data, aes(position))+geom_histogram(fill="#3300CC")
  
  if(geneInfo[,"V7"]=="+"){
  plot1 <- plot1  +
    theme_bw()+
    theme(panel.grid = element_blank(),
          panel.border = element_blank())+
    xlim(c(range[1],range[2]))+ylab("Counts Strand +")
  }else{
  plot1 <- plot1 +
    theme_bw()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank())+
    xlim(c(range[1],range[2]))+ 
    scale_y_reverse()+ylab("Counts Strand -")
} 
  plot1
  
  
}
