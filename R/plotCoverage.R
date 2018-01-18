#'@description Coverage plot for chosen part of RNA sequence
#'
#'@title plotCoverage
#'
#'@param bamDataFrame Data frame converted to R using \code{bamToR()} function.
#'@param chromosome Number of chosen chromosome.
#'@param start First position in gene/exon/three prime utr etc.
#'@param stop Last position in gene/exon/three prime utr etc.
#'@param type Type of the plot, by default \code{histogram}.
#'@param ... Optional arguments.
#'
#'@examples
#'library(directRNAExplorer)
#'data <- brmDataChromosome1[brmDataChromosome1$pos >2002000 & brmDataChromosome1$pos < 2018000,]
#'plotCoverage(data, chromosome = 1, start =2002610 , stop = 2004510)
#'
#'@importFrom ggplot2 ggplot geom_histogram theme_bw ggplotGrob aes theme element_blank labs ylab scale_y_reverse xlim geom_density
#'@importFrom grid grid.newpage grid.draw
#'@export

plotCoverage <- function(bamDataFrame, chromosome,start, stop, type="histogram",...){
  position <- bins<- NULL
  posStrand <- bamDataFrame[bamDataFrame$rname == as.character(chromosome) &
                             apply(as.data.frame(bamDataFrame$flag), 1, check_pos),
                           'pos'
                           ]
  negStrand <- bamDataFrame[bamDataFrame$rname == as.character(chromosome) &
                              apply(as.data.frame(bamDataFrame$flag), 1, check_neg),
                            'pos'
                            ]

  posStrand <- data.frame(position = posStrand)
  negStrand <- data.frame(position = negStrand)
  
  range <- c(start, stop)
  
  if(type=="histogram"){
    plot1 <- ggplot(data = posStrand , aes(position))+
      geom_histogram(col="blue", fill="blue",bins=bins)

    plot2 <- ggplot(data = negStrand , aes(position))+
      geom_histogram(col="red", fill="red", bins=bins)
  
  }
  
  if(type=="density"){
    plot1 <- ggplot(data = posStrand , aes(position))+
      geom_density(col="blue", fill="blue",stat="density")
    
    plot2 <- ggplot(data = negStrand , aes(position))+
      geom_density(col="red", fill="red",stat="density")
    
  }

    
    
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
