
#' @title calculate the Significance Score for candidate TF pairs
#' @description  function to fit time series data by nls
#' @param  correction parameters to prevent the occurrence of zero entropy due to all computational data being zero, resulting in NA (Not Available) values
#' @return  A list of gene pair, ranked by the significance score.
#' @export


two_gene_bistableswitch_circuit_simulation=function(model_para,n_replicate){
  
  
  de <- diffeqr::diffeq_setup()
  
  f <- function(u,p,t) {
    du1=p[1]*u[1]^4/(1+u[1]^4)+p[2]/(1+u[2]^4)+p[3]*u[1]
    du2=p[4]*u[2]^4/(1+u[2]^4)+p[5]/(1+u[1]^4)+p[6]*u[2]   
    
    return(c(du1,du2))
  }
  
  g <- function(u,p,t) {
    return(c(0.05,0.05)) 
  }
  
  simulation_circuit_data=function(p,n_replicate){
  
  dim1=c("gene1","gene2")
  dim2=paste("cell",1:1001,sep = "_")
  dim3=paste("replicate",1:n_replicate,sep="_")
  simulation=array(NA,dim=c(2,1001,n_replicate),dimnames = list(dim1,dim2,dim3))
  
  u0=c(0,0)
  tspan <- c(0.0,10.0)
  for ( i in 1:n_replicate){
    print(i)
    prob <- de$SDEProblem(f,g,u0,tspan,p)
    sol <- de$solve(prob,saveat=0.01)
    udf <- as.data.frame(t(sapply(sol$u,identity)))
    simulation[,,i]=t(udf)
  }
  

  return(simulation)

}

}
