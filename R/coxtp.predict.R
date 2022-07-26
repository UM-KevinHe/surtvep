
#' Title
#'
#' @param fit Model get from coxtp
#' @param baseline Baseline estimation from coxtp.baseline or the arbitary baseline entered
#' @param newdata New data in vector
#' @param strata Whether there is stratification in the dataset, default is FALSE.
#'
#' @return
#' @export
#'
#' @examples 
coxtp.predict <- function(fit,baseline,newdata=c(),strata=FALSE){
  
  
  # fit=model_result
  # baseline=baselinedata
  # newdata=data_predict
  # strata=TRUE
  
  B.spline      = as.matrix(fit$bases)
  theta         = fit$theta_list[[length(fit$theta_list)]]
  beta          = B.spline %*% t(theta)
  
  beta     <- B.spline%*%t(theta)
  var      <-  fit$VarianceMatrix
  betax   <-beta %*% newdata
  
  #No strata:
  if(!strata){
    lambda0_betax=baseline[baseline$lambda!=0,"lambda"]*exp(betax)
    outdata=data.frame(cbind(baseline[baseline$lambda!=0,"unique.time."],lambda0_betax))
    colnames(outdata)=c("unique.time.","lambda0_exp_betax") 
  } else {
    lambda0_betax=baseline[baseline$lambda!=0,"lambda"]*exp(betax)
    outdata=data.frame(cbind(baseline[baseline$lambda!=0,c("unique.time_temp.","strata")],lambda0_betax))
    colnames(outdata)=c("unique.time.","strata","lambda0_exp_betax")
    outdata$strata=as.factor(outdata$strata)
  }
  return(outdata)
  
}

