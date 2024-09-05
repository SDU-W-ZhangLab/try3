#' @title example_CNS
#' @description the example of simulation data
#' @return  A list of gene pair, ranked by the significance score.
#' @export


astrocyte_exprs=t(read.table("astrocyte.txt"))
Switch__TF=c("Hes5","Scl","Stat3","Aldh1L")
Transient_TF=c("Mash1","Zic1","Brn2","Tuj1","Olig2","Myt1L","Sox8")
SignificanceScore=get_SignificanceScore(astrocyte_exprs,Switch__TF,Transient_TF,0)
BSfate_TFs=get_singleTF_BSfate_rank(SignificanceScore)


neuron_exprs=read.table("neuron.txt")
Switch_TF=c("Mash1","Zic1","Brn2","Tuj1")
Transient_TF=c("Hes5","Scl","Stat3","Aldh1L","Olig2","Myt1L","Sox8")
SignificanceScore=get_SignificanceScore(astrocyte_exprs,Switch__TF,Transient_TF,0)
BSfate_TFs=get_singleTF_BSfate_rank(SignificanceScore)


oligodendrocyte_exprs=read.table("oligodendrocyte.txt")
Switch_TF=c("Hes5","Olig2","Myt1L","Sox8")
Transient_TF=c("Mash1","Zic1","Brn2","Tuj1","Scl","Stat3","Aldh1L")
SignificanceScore=get_SignificanceScore(astrocyte_exprs,Switch__TF,Transient_TF,0)
BSfate_TFs=get_singleTF_BSfate_rank(SignificanceScore)




layout(1)
par(mar=c(8,8,6,4))
top_predicted=cbind(1:length(SignificanceScore),SignificanceScore)

plot(top_predicted[,1],round(top_predicted[,2],3),ylim=c(-0.35,-1),xlim=c(1,28),cex = 2,col=c(rep("red",4),rep("#D3D3D3",length(top_predicted)-4)),xlab=c("Ranked candidate gene pairs"),ylab=c("Significance score"),pch=16,las=1,main="Candidate bistable switches predicted by Bsfate")
grid(lty = "dotted")
text(top_predicted[c(1:4),1]+6, top_predicted[c(1:4),2],labels = rownames(top_predicted)[c(1:4)],srt=0,cex=1.2)

layout(matrix(c(1,1, 2, 3, 4,5), ncol = 2, byrow = T), widths = c(1, 1))


par(mar=c(2,6,3,8))
plot(1,type="n",xlim=c(0,1000),ylim=c(-0.05,4.5),xlab="",ylab="EXP",frame.plot=F,las=1)
for (i in 1:length(class1)){
  points(branch_data[,class1[i]],pch=1, cex=0.5,col=brewer.pal(8, "OrRd") [8-i])

}

for(i in 1:length(class2)){
  points(branch_data[,class2[i]],pch=15, cex=0.5,col=brewer.pal(9, "Blues") [8-i+1])
}
legend(900,3,legend=c(class1,class2),col=c(brewer.pal(8, "OrRd") [c(7,6,5,4)],brewer.pal(9, "Blues") [c(8,7,6,5,4,3,2)]),lty=1,bty="n",bg="grey", ncol = 1,cex = 0.8)



bsfate_genepair=do.call(rbind,strsplit(names(entropy_score_ordered),split="_"))

par(mar=c(4,4,4,4))
for(i in 1:4){
  gene_pair=bsfate_genepair[i,]
  print(gene_pair)
  branch1_genepair=branch_data[,gene_pair]
  entropy_1=apply(branch1_genepair,1,function(x) entropy((x)/sum(x)))
  data_1=data.frame(cbind(1:length(entropy_1),entropy_1))
  colnames(data_1)=c("time","entropy")
  #model <- lm(entropy ~ time, data = data_1)
  #summary(model)
  c1=cor(data_1[,1],data_1[,2],method = "spearman",use="pairwise.complete.obs")
  print(c1)
  #plot(1,type="n",xlim=c(0,1000),ylim=c(-0.05,5),xlab="",ylab="",frame.plot=F,main="Branch1",las=1)
  #points(branch1[,gene_pair[1]],type="l",lwd=1,lty=5, col="darkred")
  #points(branch1[,gene_pair[2]],type="l",lwd=1,lty=1, col="darkgreen")
  plot(0,type="n",xlim=c(0,1000),ylim=c(0,entropy(c(0.5,0.5))),xlab="",ylab="",main=paste(c(gene_pair[1],"/",gene_pair[2]),collapse=""  ),las=1,frame.plot=F)

  points(entropy_1,pch=16, col= brewer.pal(8, "Spectral") [col_i],cex=0.8)

  text(700,0.65,paste("Score=",round(c1,2),sep=" "),cex=1.1)



}


