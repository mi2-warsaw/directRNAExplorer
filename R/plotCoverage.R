#'@description Coverage plot for chosen part of RNA sequence
#'
#'@title plotCoverage
#'
#'@param bamDataFrame Data frame converted to R using \code{bamToR()} function
#'@param chromosome Number of chosen chromosome
#'@param start First position in gene/exon/three prime utr etc.
#'@param stop Last position in gene/exon/three prime utr etc.
#'@param ... Optional arguments
#'
#'@importFrom ggplot2 ggplot geom_histogram theme_bw ggplotGrob aes theme element_blank labs ylab scale_y_reverse xlim
#'@importFrom grid grid.newpage grid.draw
#'@export

plotCoverage <- function(bamDataFrame, chromosome,start, stop, ...){
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
  plot1 <- ggplot(data = posStrand , aes(position))+
    geom_histogram(col="blue", fill="blue",bins=bins)+
    theme_bw()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank())+
    xlim(c(range[1],range[2]))
  plot1 <- plot1 + labs(title="Coverage plot")+ylab("Counts Strand +")


  plot2 <- ggplot(data = negStrand , aes(position))+
    geom_histogram(col="red", fill="red", bins=bins) +
    theme_bw()+
    theme(panel.grid = element_blank(),
      panel.border = element_blank())+
    xlim(c(range[1],range[2]))

  plot2 <- plot2 + scale_y_reverse()+ylab("Counts Strand -")


  grid.newpage()
  grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))
  }
