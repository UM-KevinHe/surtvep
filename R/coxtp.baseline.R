coxtp.baseline <- function(fit, delta,z,time){
  

  
  model1  = fit
  
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
  
  return(baselinedata)
  
}