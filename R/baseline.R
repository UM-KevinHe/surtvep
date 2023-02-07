#' calculating baseline hazard and baseline cumulative hazard using the result from a `coxtv` or `coxtp` object
#'
#' The baseline estimation is the baseline hazard at each observed failure time when 
#' holding all the covariates to be zero.
#'
#' @param fit model from `coxtv` or `coxtp`.
#' 
#' @export
#' 
#' @return A list with three components:
#' \item{time}{the unique observed failure time.} 
#' \item{hazard}{the baseline hazard corresponding to each unique failure time point.}
#' \item{cumulHaz}{the cumulative baseline hazard corresponding to each unique failure time point.}
#'
#' @examples
#' \dontrun{
#' data(ExampleData)
#' z <- ExampleData$x
#' time  <- ExampleData$time
#' event <- ExampleData$event
#' 
#' fit   <- coxtv(event = event, z = z, time = time)
#' base.est <- baseline(fit)
#' }
#' 
#' 
baseline <-  function(fit, ...){
  
    
  if (missing(fit)) stop ("Argument fit is required!")
  if (class(fit)!="coxtp" & class(fit)!="coxtv") stop("Object fit is not of class 'coxtv' or 'coxtp'!")
  
  
  
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
  k=ncol(fit$bases)
  result1 <- Lambda_estimate_ties2(knot = k, delta = event,
                                   z = as.matrix(z), b_spline = as.matrix(B.spline),
                                   theta = theta_IC, tieseq = tieseq)
  lambda   <- result1$lambda
  time2 <- unique(time)
  Lambda <- cumsum(lambda)
  
  baselinedata <- data.frame(unique(time),lambda,Lambda)
  # colnames(baselinedata) <- c("time", "hazard", "Lambda")
  
  res <- list("time" = unique(time),
              "hazard" = as.numeric(lambda),
              "cumulHaz" = Lambda)
  
  class(res) <- "baseline"
  
  return(res)
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






#' plotting baseline
#' 
#' Plotting baseline from a fitted `baseline` object.
#'
#' @param fit model get from `baseline` function.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon theme_bw theme element_text element_blank margin labs ggtitle
#'
#' @exportS3Method plot baseline
#'
#' @examples
#' \dontrun{
#' data(ExampleData)
#' z <- ExampleData$x
#' time  <- ExampleData$time
#' event <- ExampleData$event
#' 
#' fit   <- coxtv(event = event, z = z, time = time)
#' base.est <- baseline(fit)
#' }
plot.baseline <- function(fit, xlab, ylab, xlim, ylim,
                                title,...){
  
  
  if (missing(xlab)) xlab <- "time"
  if (missing(ylab)) ylab <- "cumulative hazard"
  
  if (missing(fit)) stop ("Argument fit is required!")
  if (class(fit)!="baseline") stop("Object fit is not of class 'baseline'!")
  
  missingxlim <- missing(xlim); missingylim <- missing(ylim); 
  missingtitle <- missing(title);
  
  plot_data <- data.frame(fit$time, fit$cumulHaz)
  colnames(plot_data) <- c("time", "cumulHaz")
  
  plt <- ggplot(data = plot_data, aes(x = time, y = cumulHaz)) +
    geom_line( size = 0.9)+
    labs(x="time", y = "cumulative hazard")+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5))

  if (missingxlim) {
    plt <- plt + scale_x_continuous(name=xlab)
  } else {
    if (!is.numeric(xlim)) stop("Invalid xlim!")
    plt <- plt + scale_x_continuous(name=xlab, limits=xlim)
  }
  if (missingylim) {
    plt <- plt +
      scale_y_continuous(name=ylab)
  } else {
    if (!is.numeric(ylim)) stop("Invalid ylim!")
    plt <- plt +
      scale_y_continuous(name=ylab, limits=ylim)
  }
  if (!missingtitle) plt <- plt + ggtitle(title)
  
  return(plt)
  
}
