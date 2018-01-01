#'@description Coverage plot for chosen part of RNA sequence
#'
#'@title plotCoverage
#'
#'@param bamDataFrame data frame converted to R using \code{bamToR()} function
#'@param chromosome number of chosen chromosome
#'@param ... optional arguments
#'
#'@importFrom ggplot2 ggplot geom_histogram theme_bw ggplotGrob aes theme element_blank labs ylab scale_y_reverse
#'@importFrom grid grid.newpage grid.draw
#'@export

plotCoverage <- function(bamDataFrame, chromosome,...){
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

  plot1 <- ggplot(data = posStrand , aes(position))+
    geom_histogram(col="blue", fill="blue",bins=bins)+
    theme_bw()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank())
  plot1 <- plot1 + labs(title="Coverage plot")+ylab("Counts Strand +")


  plot2 <- ggplot(data = negStrand , aes(position))+
    geom_histogram(col="red", fill="red", bins=bins) +
    theme_bw()+
    theme(panel.grid = element_blank(),
      panel.border = element_blank())

  plot2 <- plot2 + scale_y_reverse()+ylab("Counts Strand -")


  grid.newpage()
  grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))
  }
