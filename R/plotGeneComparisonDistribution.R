#'@title plotGeneComparisonDistribution
#'
#'@description Plot of distribution of two samples on selected gene
#'
#'@param gene Character with gene name.
#'@param data1 Data frame converted to R using \code{bamToR()} function, corresponding to the first sample.
#'@param data2 Data frame converted to R using \code{bamToR()} function, corresponding to the second sample.
#'@param geneData Data frame with positions of all genes and their names, by default we use a \code{TAIR10_genes}.
#'@param range How many nucleotide before \code{start} and after \code{stop} we include to genes.
#'@param genePart The part of gene we want to visualize.
#'@param adjust Adjust parameter for the density kernel
#'@param ... Optional arguments.
#'
#'@importFrom dplyr filter
#'@importFrom tidyr separate
#'
#' @export


plotGeneComparisonDistribution <- function(gene, data1, data2, geneData=SequencingExplainer::TAIR10_genes, range=0, genePart="gene", adjust = 1,...){
  pos<-id<-V4<-V5<-V3 <- NULL
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
  
  bam1<- data1@listData
  bam1 <- bam1[c(1:11)]
  bam1 <- as.data.frame(bam1)
  bamDataFrameGene1 <- filter(bam1, pos>=startPosition & pos <=endPosition)
  bamDataFrameGene1 <- bamDataFrameGene1[which(bamDataFrameGene1$rname==chromosome),]
  
  bam2<- data2@listData
  bam2 <- bam2[c(1:11)]
  bam2 <- as.data.frame(bam2)
  bamDataFrameGene2 <- filter(bam2, pos>=startPosition & pos <=endPosition)
  bamDataFrameGene2 <- bamDataFrameGene2[which(bamDataFrameGene2$rname==chromosome),]
  
  ####
  
  posStrand1 <- bamDataFrameGene1[bamDataFrameGene1$rname == as.character(chromosome) &
                              apply(as.data.frame(bamDataFrameGene1$flag), 1, check_pos),
                            'pos'
                            ]
  negStrand1 <- bamDataFrameGene1[bamDataFrameGene1$rname == as.character(chromosome) &
                              apply(as.data.frame(bamDataFrameGene1$flag), 1, check_neg),
                            'pos'
                            ]
  
  posStrand1 <- data.frame(position = posStrand1)
  negStrand1 <- data.frame(position = negStrand1)
  
  posStrand2 <- bamDataFrameGene2[bamDataFrameGene2$rname == as.character(chromosome) &
                              apply(as.data.frame(bamDataFrameGene2$flag), 1, check_pos),
                            'pos'
                            ]
  negStrand2 <- bamDataFrameGene2[bamDataFrameGene2$rname == as.character(chromosome) &
                              apply(as.data.frame(bamDataFrameGene2$flag), 1, check_neg),
                            'pos'
                            ]
  
  posStrand2 <- data.frame(position = posStrand2)
  negStrand2 <- data.frame(position = negStrand2)
  
  range <- c(startPosition, endPosition)
  

  
    plot1 <- ggplot(data = posStrand1 , aes(position))+
      geom_density(col="#CC0033",stat="density", adjust = adjust)+
      geom_density(data = posStrand2,col="#3300CC", stat="density", adjust = adjust)
    
    plot2 <- ggplot(data = negStrand1 , aes(position))+
      geom_density(col="#CC3300",stat="density", adjust = adjust)+
      geom_density(data = negStrand2,col="#0033CC", stat="density", adjust = adjust)
    

  
  
  plot1 <- plot1  +
    theme_bw()+
    theme(panel.grid = element_blank(),
          panel.border = element_blank())+
    xlim(c(range[1],range[2]))+
    labs(title="Coverage plot")+ylab("Counts Strand +")
  
  plot2 <- plot2 +
    theme_bw()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank())+
    xlim(c(range[1],range[2]))+ 
    scale_y_reverse()+ylab("Counts Strand -")
  
  
  grid.newpage()
  grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))
  invisible(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))
  
  
  
}