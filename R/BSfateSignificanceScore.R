#' @title calculate the Significance Score for candidate TF pairs
#' @description  function to fit time series data by nls
#' @param  correction parameters to prevent the occurrence of zero entropy due to all computational data being zero, resulting in NA (Not Available) values
#' @return  A list of gene pair, ranked by the significance score.
#' @export


get_SignificanceScore=function(scExp,tf_1,tf_2,entropy_correction_para){

  combinations_genepairs <- data.frame(expand.grid(unique(tf_1),unique(tf_2)))

  SignificanceScore=c()
  for( i in 1:nrow(combinations_genepairs)){

    gene1=as.character(combinations_genepairs[i,1])
    gene2=as.character(combinations_genepairs[i,2])

    exp_genepair=t(scExp[c(gene1,gene2),])

    entropy_pairwise=apply(exp_genepair,1,function(x) entropy((x+entropy_correction_para)/sum(x+entropy_correction_para)))
    SC=cor(1:length(entropy_pairwise),entropy_pairwise,method = "spearman",use = "pairwise.complete.obs")

    SignificanceScore=c(SignificanceScore,SC)
  }
  names(SignificanceScore)=apply(combinations_genepairs,1,function(x) paste(x,collapse = "_"))
  SignificanceScore=SignificanceScore[order(SignificanceScore)]

  SignificanceScore=SignificanceScore[!is.na(SignificanceScore)]
  return(SignificanceScore)

}


