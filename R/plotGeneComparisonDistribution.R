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
#'@param adjust A multiplicate bandwidth adjustment.
#'@param stat A type of plot.
#'@param geneName A name of presented gene.
#'@param ... Optional arguments.
#'
#'@importFrom dplyr filter
#'@importFrom tidyr separate
#'@importFrom ggplot2 ggplot scale_fill_manual theme_bw ggplotGrob aes theme element_blank labs ylab scale_y_reverse xlim stat_ecdf scale_color_manual geom_hline element_text
#'@importFrom grid grid.newpage grid.draw
#' @export


plotGeneComparisonDistribution <- function(gene, data1, data2, geneData=directRNAExplorer::TAIR10_genes, range=0, genePart="gene", adjust = 1, stat="density", geneName = NULL,...){
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
  
 
  bamDataFrameGene1 <- data1[which(data1$pos>=startPosition & data1$pos <=endPosition),]
#bamDataFrameGene1 <- subset(data1, pos>=startPosition & pos <=endPosition)
  bamDataFrameGene1 <- bamDataFrameGene1[which(bamDataFrameGene1$rname==chromosome),]
  
  
  bamDataFrameGene2 <- data2[which(data2$pos>=startPosition & data2$pos <=endPosition),]
  bamDataFrameGene2 <- bamDataFrameGene2[which(bamDataFrameGene2$rname==chromosome),]
  
  
  position <- c(bamDataFrameGene1$pos, bamDataFrameGene2$pos)
  strand <- c(as.character(bamDataFrameGene1$strand), as.character(bamDataFrameGene2$strand))
  condition <- c(rep("type1", length(bamDataFrameGene1$pos)), rep("type2", length(bamDataFrameGene2$pos)))
  
  data <- data.frame(position, strand, condition)
  
  posStrand <- dplyr::filter(data, strand=="+")
  negStrand <- dplyr::filter(data, strand == "-")
  range <- c(startPosition, endPosition)
  
  if(stat=="density"){
  
    plot1 <- ggplot(data = posStrand , aes(position, col=condition))+
      geom_density(stat="density", adjust = adjust, show_guide=TRUE)
      #geom_density(data = posStrand2,col="#3300CC", stat="density", adjust = adjust, show.legend=TRUE)
    
    plot2 <- ggplot(data = negStrand , aes(position, col=condition))+
      geom_density(col="#CC3300",stat="density", adjust = adjust, show.legend=TRUE)
      #geom_density(data = negStrand2,col="#0033CC", stat="density", adjust = adjust, show.legend=TRUE)
  }
  
  if(stat=="ecdf"){
    
    plot1 <- ggplot(posStrand, aes(position, col=condition))+
      stat_ecdf(show.legend=TRUE)+scale_color_manual(values=c("#CC0033", "#3300CC"))
    plot2 <- ggplot(negStrand, aes(position, col=condition))+stat_ecdf(show.legend=TRUE)+scale_color_manual(values=c("#CC3300", "#0033CC"))#scale_color_manual(values = c("#3300CC", "#0033CC"))
    
  }

  
  
  plot1 <- plot1  +
    theme_bw()+
    theme(panel.grid = element_blank(),
          panel.border = element_blank())+
    xlim(c(range[1],range[2]))+ylab("Counts Strand +")
  
  plot2 <- plot2 +
    theme_bw()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank())+
    xlim(c(range[1],range[2]))+ 
    scale_y_reverse()+ylab("Counts Strand -")
  
  if(!is.null(geneName)){
    plot1 <- plot1 + theme(plot.title = element_text(paste0("Gene ", geneName)))
    plot2 <- plot2 + theme(plot.title = element_text(paste0("Gene ", geneName)))
  }
  
  
  if(nrow(plot2$data)==0){
  grid.newpage()
  grid.draw(rbind(ggplotGrob(plot1), size = "last"))
  invisible(rbind(ggplotGrob(plot1), size = "last"))
  
  }else if(nrow(plot1$data)==0){
    grid.newpage()
    grid.draw(rbind(ggplotGrob(plot2), size = "last"))
    invisible(rbind(ggplotGrob(plot2), size = "last"))
    
  }else{
    grid.newpage()
    grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))
    invisible(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))
    
  }
  
}