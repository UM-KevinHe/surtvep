#' Calculating Baseline hazard using the result from coxtp()
#'
#' @param fit Model get from coxtp
#' @param delta event vector, should be a vector containing 0 or 1
#' @param z Covariate matrix
#' @param time Time vector, should be a vector with non-negative numeric value
#' @param strata stratification group defined in the data. If there exist stratification group, please enter as vector.
#'
#' @return
#' @export
#'
#' @examples
coxtp.baseline <- function(fit, delta,z,time,strata=c()){
  

  
  model1  = fit
  # 
  # fit=model1$model_result
  # delta=event
  # z=data
  # time=round(time,2)
  # strata=strata
  # 
  # model1=model1$model_result
  
  if(length(strata)==0){
    #Calculated unique time and ties for the baseline calculation:
    unique_time   <- unique(time)
    
    tieseq <- NULL
    index  <- NULL
    for (i in 1:length(unique(time))) {
      tieseq[i] <- length(which(time==unique(time)[i]))
      index[[i]]  <- (which(time==unique(time)[i]))
    }
    
    #call Rcpp
    theta_IC <- model1$theta_list[[length(model1$theta_list)]]
    B.spline=bSpline(unique_time,knots=model1$knots,intercept=TRUE,degree=3)
    k=ncol(model1$bases)
    result1 <- Lambda_estimate_ties2(knot = k, delta = delta,
                                     z = z, b_spline = B.spline,
                                     theta = theta_IC, tieseq = tieseq)
    lambda   <- result1$lambda
    time2 <- unique(time)
    Lambda <- cumsum(lambda)
    
    baselinedata <- data.frame(unique(time),lambda,Lambda)
  } else {
    flevel=unique(strata)
    for(f in flevel){
      time_temp=time[strata==f]
      event_temp=delta[strata==f]
      data_temp=z[strata==f,]

      #Calculated unique time and ties for the baseline calculation:
      unique_time   <- unique(time_temp)
      tieseq <- NULL
      index  <- NULL
      for (i in 1:length(unique(time_temp))) {
        tieseq[i] <- length(which(time_temp==unique(time_temp)[i]))
        index[[i]]  <- (which(time_temp==unique(time_temp)[i]))
      }

      #call Rcpp
      theta_IC <- model1$theta_list[[length(model1$theta_list)]]
      B.spline=bSpline(unique_time,knots=model1$knots,intercept=TRUE,degree=3)
      k=ncol(model1$bases)
      result1 <- Lambda_estimate_ties2(knot = k, delta = event_temp,
                                       z = data_temp, b_spline = B.spline,
                                       theta = theta_IC, tieseq = tieseq)
      lambda   <- result1$lambda
      time2 <- unique(time_temp)
      Lambda <- cumsum(lambda)

      temp_data <- data.frame(unique(time_temp),lambda,Lambda)
      temp_data$strata=f
      if(f!=flevel[1]){
        baselinedata=rbind(baselinedata,temp_data)
      } else {
        baselinedata=temp_data
      }

    }
  }
  
  return(baselinedata)
  
}
