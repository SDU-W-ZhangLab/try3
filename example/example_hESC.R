#' @title example_hESC
#' @description the example of hESC
#' @return  A list of gene pair, ranked by the significance score.
#' @export


library(stats)
library(entropy)
library("RColorBrewer")

load("scExp_hESC.RData")
load("pesudo_hESC.RData")
TF_Human=read.csv("TF_human.txt")[,1]


scExp_hESC_TF=scExp_hESC[intersect(rownames(scExp_hESC),TF_Human),order(pesudo_hESC[,1],decreasing = F)]

Switch_test_hESC=Switch_nls(scExp_hESC_TF)
Tansient_test_hESC=Tansient_nls(scExp_hESC_TF)

Pro_screen_TFs=Screen_TF(Switch_test_hESC, Tansient_test_hESC,top_n)
SignificanceScore=get_SignificanceScore(scExp_hESC_TF,Hs_switch_n,Hs_transient_n,0.1)
BSfate_TFs=get_singleTF_BSfate_rank(SignificanceScore)









