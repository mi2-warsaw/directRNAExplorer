#'plotPoints 
#'
#'@param gene Character with gene name.
#'@param data1 Data frame converted to R using \code{bamToR()} function.
#'@param data2 Data frame converted to R using \code{bamToR()} function.
#'@param data3 Data frame converted to R using \code{bamToR()} function.
#'@param geneData Data frame with positions of all genes and their names, by default we use a \code{Araport11}.
#'@param range_st How many nucleotides before \code{start} include to genes.
#'@param range_end How many nucleotides after \code{stop} we include to genes.
#'
#'
#' @importFrom dplyr filter group_by summarise 
#' @importFrom ggplot2 ggplot geom_point facet_grid aes element_text
#' 
#' @export

plotPoints <- function(gene, data1, data2, data3, geneData = directRNAExplorer::Araport11, range_st=0 , range_end=0){
  rname<-V3<-id<-pos<-count<-type<-position<-number<-NULL
  
  dic_gene <- filter(geneData, V3=="gene")
  gene_data <- filter(dic_gene, id == gene)
  
  start <- gene_data$V4
  stop <- gene_data$V5
  chromosome <- gene_data$V1
  chromosome <- substr(chromosome, 4,4)
  strand <- factor(gene_data$V7)
  
  if(strand=="+"){
    data1 <- data1[which(data1$strand=="-"),]
    data2 <- data2[which(data2$strand=="-"),]
    data3 <- data3[which(data3$strand=="-"),]
  }else{
    data1 <- data1[which(data1$strand=="+"),]
    data2 <- data2[which(data2$strand=="+"),]
    data3 <- data3[which(data3$strand=="+"),]
  }
  
  data1 <- data1[which(data1$rname==chromosome),]
  data2 <- data2[which(data2$rname==chromosome),]
  data3 <- data3[which(data3$rname==chromosome),]
  
  pos1 <- data1$pos
  type1 <-  rep(factor(1), length(data1$pos))
  
  data1 <- data.frame(pos=pos1, type=type1)
  
  pos2 <- data2$pos
  type2 <-  rep(factor(2), length(data2$pos))
  
  data2 <- data.frame(pos=pos2, type=type2)
  
  
  pos3 <- data3$pos
  type3 <-  rep(factor(3), length(data3$pos))
  
  data3 <- data.frame(pos=pos3, type=type3)
  
  if(range_st >0 || range_end >0){
  data1_plot <- filter(data1, pos>=(start-range_st) & pos<=(stop+range_end))
  data2_plot <- filter(data2, pos>=(start-range_st) & pos<=(stop+range_end))
  data3_plot <- filter(data3, pos>=(start-range_st) & pos<=(stop+range_end))
  }else{
    data1_plot <- filter(data1, pos>=start & pos<=stop)
    data2_plot <- filter(data2, pos>=start & pos<=stop)
    data3_plot <- filter(data3, pos>=start & pos<=stop)
  }
  
  
  dataAll <- rbind(data1_plot, data2_plot, data3_plot)
  dataAll <- as.data.frame(dataAll)
  dataAll$type <- factor(dataAll$type)
  colnames(dataAll)[1] <- "position"
  
  dataAll$number <- rep(1, nrow(dataAll))
  groupedData <- group_by(dataAll, position, type)
  groupedData <- summarise(groupedData, count = sum(number))
  

  plot <- ggplot(groupedData, aes(position,count, col=type))+
   geom_point()+
  facet_grid(type~.)
  
  return(plot)
  
}