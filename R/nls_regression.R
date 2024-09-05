#' @title screen switch-like and transient TFs
#' @description  function to fit time series data by nls
#' @param lineage1 refer to time series of committed data
#' @return  G1 & G2
#' @export

Switch_nls<-function(scExp){

  nls_fit<-function(Data){

    data <- list(
      y = as.vector(Data),
      x = (1:length(Data)/length(Data))
    )

    nls_model <- try(nls(y ~ theta1 / (1 + exp(-10*theta2 * (x - theta3))),data = data, start = list( theta1 = max(Data), theta2 = 1, theta3 = 0.5)),silent = T)

    if(!'try-error' %in% class(nls_model))
    {

      params <- coef(nls_model)
      fitted_values=predict(nls_model)
      RSS <- sum((data[["y"]] - fitted_values)^2)
      TSS <- sum((data[["y"]] - mean(data[["y"]]))^2)
      rsquared <- 1 - RSS / TSS

    }else{
      params=rep(NA,3)
      rsquared=NA
    }
    model_summary=c(rsquared,params)
    names(model_summary)=c("R-squared","Asym", "scal","xmid")
    return(model_summary)
  }

  switch_test=t(apply(scExp,1,nls_fit))
  switch_test=as.data.frame(na.omit(switch_test))

  return(switch_test)}



Tansient_nls<-function(scExp){

  nls_fit<-function(Data){

    data <- list(
      y = as.vector(Data),
      x = (1:length(Data)/length(Data))
    )

    nls_model <- try(nls(y ~ theta1 *exp(-10*theta2 * (x - theta3)^2),data = data, start = list(theta1 = max(Data), theta2 = 5, theta3 = 0.5)),silent = T)



    if(!'try-error' %in% class(nls_model))            # 判断当前循环的try语句中的表达式是否运行正确
    {

      params <- coef(nls_model)
      fitted_values=predict(nls_model)

      RSS <- sum((data[["y"]] - fitted_values)^2)
      TSS <- sum((data[["y"]] - mean(data[["y"]]))^2)
      rsquared <- 1 - RSS / TSS

    }else{
      params=rep(NA,3)
      rsquared=NA
    }
    model_summary=c(rsquared,params)
    names(model_summary)=c("R-squared","Asym", "scal","xmid")
    return(model_summary)
  }


  transient_test=t(apply(scExp,1,nls_fit))
  transient_test=as.data.frame(na.omit(transient_test))

  return(transient_test)}




