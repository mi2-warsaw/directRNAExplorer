#'@description Converts .bam file into a R list
#'
#'@title bamToRs
#'
#'@param filename string contatinig name of file with bam data. Remember that you need to set path to folder with selected file.
#'
#'
#'@importFrom Rsamtools scanBam
#'
#'@examples
#'\dontrun{
#'filename <- "P13.2016-11-17T16_40_36.1234062861011.fc1.ch1.sorted.bam"
#'bamData <- bamToR(filename)
#'}
#'@export


bamToR <- function(filename){
  if(sub('.*\\.', '', filename)!="bam"){
    stop("Please specify the .bam file")
  }

  bamFile <- scanBam(filename)

  .unlist <- function (x){
    x1 <- x[[1L]]
    if (is.factor(x1)){
      structure(unlist(x), class = "factor", levels = levels(x1))
    } else {
      do.call(c, x)
    }
  }

  bamField <- names(bamFile[[1]])

  list <- lapply(bamField, function(y) .unlist(lapply(bamFile, "[[", y)))

  bamDataFrame <- do.call("DataFrame", list)
  names(bamDataFrame) <- bamField

  return(bamDataFrame)
}
