
#' @title screen switch-like TFs
#' @description  function to fit time series data by nls
#' @param lineage1 refer to time series of committed data
#' @return  G1 & G2
#' @export


Screen_TF<-function(Switch_test, Tansient_test,top_n){
  
  Switch_test=Switch_test[Switch_test[,"scal"]>0.6&Switch_test[,"scal"]<15,]
  Switch_test=Switch_test[order(Switch_test[,"R-squared"],decreasing = T),]
  
  Switch_test_top_n=Switch_test[1:min(top_n,nrow(Switch_test)),]
  
  Tansient_test=Tansient_test[which(Tansient_test[,"xmid"]<= 0.7&Tansient_test[,"xmid"]>0.3&Tansient_test[,"scal"]>0.6&Tansient_test[,"scal"]<10&Tansient_test[,"Asym"]>0.2),]
  Tansient_test=Tansient_test[order(Tansient_test[,"R-squared"],decreasing = T),]
  Tansient_test_top_n=Tansient_test[1:min(top_n,nrow(Tansient_test)),]
  
  TF_candidate=cbind(rownames(Switch_test_top_n),rownames(Tansient_test_top_n))
  
  return(TF_candidate)}


