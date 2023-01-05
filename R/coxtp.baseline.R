#' calculating baseline hazard using the result from a `coxtv` object
#'
#' The baseline estimation is the baseline hazard at each observation time when holding all the covariates equals to zero.
#'
#' @param fit model from `coxtp`
#' 
#' @export
#' 
#' @return a list with three components
#' \item{time}{the unique observation time} 
#' \item{hazard}{the baseline hazard corresponding to each unqiue time point}
#' \item{cumulHaz}{the cumulative baseline hazard corresponding to each unqiue time point}
#'
#' @examples
#' \dontrun{
#' data(ExampleData)
#' z <- ExampleData$x
#' time  <- ExampleData$time
#' event <- ExampleData$event
#' 
#' lambda  = c(0,1)
#' fit   <- coxtp(event = event, z = z, time = time, lambda=lambda)
#' base.est <- coxtp(fit)
#' }
#' 
#' 
coxtp.baseline <-  function(fit, ...){
  
    
  if (missing(fit)) stop ("Argument fit is required!")
  if (class(fit)!="coxtp") stop("Object fit is not of class 'coxtp'!")
  
  model1  = fit
  data       <- attr(fit, "data")
  event  <- data$event
  strata <- data$strata
  time   <- data$time
  z      <- subset(data, select = -c(event, strata, time))
  
  
  
  
  # if(length(strata)==0){
    #Calculated unique time and ties for the baseline calculation:
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
    k=ncol(model1$bases)
    result1 <- Lambda_estimate_ties2(knot = k, delta = event,
                                     z = as.matrix(z), b_spline = as.matrix(B.spline),
                                     theta = theta_IC, tieseq = tieseq)
    lambda   <- result1$lambda
    time2 <- unique(time)
    Lambda <- cumsum(lambda)
    
    baselinedata <- data.frame(unique(time),lambda,Lambda)
    # colnames(baselinedata) <- c("time", "hazard", "Lambda")
    return(list("time" = unique(time),
                "hazard" = as.numeric(lambda),
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






