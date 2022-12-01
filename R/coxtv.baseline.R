#' calculating Baseline hazard using the result from a `coxtv` object
#'
#' @param fit Model get from coxtp
#' @param delta event vector, should be a vector containing 0 or 1
#' @param z Covariate matrix
#' @param time Time vector, should be a vector with non-negative numeric value
#' @param strata stratification group defined in the data. If there exist stratification group, please enter as vector.
#'
#' @export
#' 
#' @return a list with components
#' \item{time}{the unique event time points} 
#' \item{hazard}{the baseline hazard corresponding to each unqiue time point}
#' \item{cumulHaz}{the cumulative baseline hazard corresponding to each unqiue time point}
#' 
#' @examples
#'data("ExampleDataBinary")
#'z <- ExampleDataBinary$x
#'time <- ExampleDataBinary$time
#'event <- ExampleDataBinary$event
#'fit <- coxtv(event = event, z = z, time = time)
#'base.est = coxtv.baseline(fit, event, z, time)


coxtv.baseline <- function(fit, event, z, time, strata= NULL, ...){
  if (missing(fit)) stop ("Argument fit is required!")
  if (class(fit)!="coxtv") stop("Object fit is not of class 'coxtv'!")
  
  # if(length(strata)==0){
    unique_time   <- unique(time)
    
    tieseq <- NULL
    index  <- NULL
    for (i in 1:length(unique(time))) {
      tieseq[i] <- length(which(time==unique(time)[i]))
      index[[i]]  <- (which(time==unique(time)[i]))
    }
    
    #call Rcpp
    theta_IC <- fit$ctrl.pts
    B.spline=splines::bs(unique_time,knots=fit$internal.knots,intercept=TRUE,degree=3)
    k = ncol(fit$bases)
    result1 <- Lambda_estimate_ties2(knot = k, delta = event,
                                     z = z, b_spline = B.spline,
                                     theta = theta_IC, tieseq = tieseq)
    lambda   <- result1$lambda
    time2 <- unique(time)
    Lambda <- cumsum(lambda)
    
    baselinedata <- data.frame(unique(time),lambda,Lambda)
    # colnames(baselinedata) <- c("time", "hazard", "Lambda")
    return(list("time" = unique(time),
                "hazard" = lambda,
                "cumulHaz" = Lambda))
    
  # } else {
  #   flevel=unique(strata)
  #   for(f in flevel){
  #     time_temp=time[strata==f]
  #     event_temp=delta[strata==f]
  #     data_temp=z[strata==f,]
  #     
  #     #Calculated unique time and ties for the baseline calculation:
  #     unique_time   <- unique(time_temp)
  #     tieseq <- NULL
  #     index  <- NULL
  #     for (i in 1:length(unique(time_temp))) {
  #       tieseq[i] <- length(which(time_temp==unique(time_temp)[i]))
  #       index[[i]]  <- (which(time_temp==unique(time_temp)[i]))
  #     }
  #     
  #     #call Rcpp
  #     theta_IC <- model1$theta_list[[length(model1$theta_list)]]
  #     B.spline=bSpline(unique_time,knots=model1$knots,intercept=TRUE,degree=3)
  #     k=ncol(model1$bases)
  #     result1 <- Lambda_estimate_ties2(knot = k, delta = event_temp,
  #                                      z = data_temp, b_spline = B.spline,
  #                                      theta = theta_IC, tieseq = tieseq)
  #     lambda   <- result1$lambda
  #     time2 <- unique(time_temp)
  #     Lambda <- cumsum(lambda)
  #     
  #     temp_data <- data.frame(unique(time_temp),lambda,Lambda)
  #     temp_data$strata=f
  #     if(f!=flevel[1]){
  #       baselinedata=rbind(baselinedata,temp_data)
  #     } else {
  #       baselinedata=temp_data
  #     }
  #     
  #   }
  # }
  
  # return(baselinedata)
  
}
