#' @title example_mESC
#' @description the example of mESC
#' @return  A list of gene pair, ranked by the significance score.
#' @export

library(stats)
library(entropy)
library("RColorBrewer")

load("scExp_mESC.RData")
load("pesudo_mESC.RData")
TF_Mouse=read.csv("TF_mouse.txt")[,1]


scExp_mESC_TF=scExp_mESC[intersect(rownames(scExp_mESC),TF_Mouse),order(pesudo_mESC[,1],decreasing = F)]
Switch_test_mESC=Switch_nls(scExp_mESC_TF)
Tansient_test_mESC=Tansient_nls(scExp_mESC_TF)
top_n=20
Pro_screen_TFs=Screen_TF(Switch_test_mESC, Tansient_test_mESC,top_n)
SignificanceScore=get_SignificanceScore(scExp_mESC_TF,Pro_screen_TFs[,1],Pro_screen_TFs[,2],0.01)
BSfate_TFs=get_singleTF_BSfate_rank(SignificanceScore)






