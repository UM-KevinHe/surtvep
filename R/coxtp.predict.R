
#' make predictions from a "coxtv" or "coxtp" object.
#' 
#' @description 
#' Similar to other predict methods, this functions predicts fitted values, logits, coefficients and more from a fitted ""coxtv" or "coxtp"object.
#' 
#' @param fit Model get from coxtp
#' @param baseline Baseline estimation from coxtp.baseline or the arbitary baseline entered
#' @param newdata New data in vector
#' @param strata Whether there is stratification in the dataset, default is FALSE.
#' @param out_seq Output time sequence. Default is out_seq=NULL. When default, will output the time sequence according to input fit
#'
#' @return
#' @export
#'
#' @examples 
coxtp.predict <- function(fit,baseline,newdata=c(),strata=FALSE,out_seq=NULL){
  
  
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
  if(is.null(out_seq)){
    return(outdata)
  } else {
    for(i in 1:length(out_seq)){
      out_time=out_seq[i]
      sort_time=sort(c(outdata$unique.time.,out_time))
      out_position=which(sort_time==out_time)
      #Decide points to calculate line function
      if(out_position==1){
        pre_value=outdata[out_position,]
        next_value=outdata[out_position+1,]
      } else if(out_position>nrow(outdata)) {
        pre_value=outdata[out_position-2,]
        next_value=outdata[out_position-1,]
        
      } else {
        pre_value=outdata[out_position-1,]
        next_value=outdata[out_position,]
      }
      #Calculate line function
      slope=(next_value[1,2]-pre_value[1,2])/(next_value[1,1]-pre_value[1,1])
      b=next_value[1,2]-slope*next_value[1,1]
      predict_value=slope*out_time+b
      #merge data
      if(i==1){
        outdata_seq=c(out_time,predict_value)
      } else {
        outdata_seq=rbind(outdata_seq,c(out_time,predict_value))
      }
    }
    #rename variable
    colnames(outdata_seq)=colnames(outdata)
    rownames(outdata_seq)=NULL
    outdata_seq=as.data.frame(outdata_seq)
    return(outdata_seq)
  }
  
}

