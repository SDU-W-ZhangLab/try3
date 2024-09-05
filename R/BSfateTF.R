
#' @title get_single_BSfate_rank 
#' @description  rank single tf according to the rank of candidate TF pairs
#' @return  A list of TF, ranked by the significance score.
#' @export


get_singleTF_BSfate_rank=function(SignificanceScore){
  
  interleave_vectors <- function(a, b) {
    max_length <- max(length(a), length(b))
    
    if (length(a) < max_length) {
      a <- c(a, rep(NA, max_length - length(a)))
    }
    if (length(b) < max_length) {
      b <- c(b, rep(NA, max_length - length(b)))
    }
    
    c <- vector("numeric", 2 * max_length)
    c[seq(1, length(c), 2)] <- a
    c[seq(2, length(c), 2)] <- b
    
    c <- c[!is.na(c)]
    
    return(c)
  }
  
  BSfate_Split=do.call(rbind, strsplit(names(SignificanceScore),,split="_"))
  BSfate_TF_1=unique(BSfate_Split[,1])
  BSfate_TF_2=unique(BSfate_Split[,2])
  BSfate_TFs <- unique(interleave_vectors(BSfate_TF_1, BSfate_TF_2))
  return(BSfate_TFs)
  
}
