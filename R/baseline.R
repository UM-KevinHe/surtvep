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
#' data(ExampleData)
#' z <- ExampleData$x
#' time  <- ExampleData$time
#' event <- ExampleData$event
#' 
#' fit   <- coxtv(event = event, z = z, time = time)
#' base.est <- baseline(fit)
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
  
  event  <- event[time_order <- order(time)]
  z      <- z[time_order,]
  strata <- strata[time_order]
  time   <- time[time_order]
  
  base2 <- attr(fit, "basehazard")
  base2 <- as.numeric(base2[[1]])
  unique_time_event = unique(time[event==1])
  unique_time <- unique(time)
  
  lambda <- rep(0, length(unique_time))
  for (i in 1:length(unique_time_event)) {
    time_tmp <- unique_time_event[i]
    lambda[unique_time==time_tmp] <- base2[i]
  }

  Lambda <- cumsum(lambda)
  
  
  baselinedata <- data.frame(unique_time,lambda,Lambda)
  
  res <- list("time" = unique(time),
              "hazard" = as.numeric(lambda),
              "cumulHaz" = Lambda)
  
  class(res) <- "baseline"
  
  return(res)
  
  # if (missing(fit)) stop ("Argument fit is required!")
  # if (class(fit)!="coxtp" & class(fit)!="coxtv") stop("Object fit is not of class 'coxtv' or 'coxtp'!")
  # 
  # data       <- attr(fit, "data")
  # event  <- data$event
  # strata <- data$strata
  # time   <- data$time
  # z      <- subset(data, select = -c(event, strata, time))
  # 
  # event  <- event[time_order <- order(time)]
  # z      <- z[time_order,]
  # strata <- strata[time_order]
  # time   <- time[time_order]
  # stratum <- if (length(strata) == 0) rep(1, length(time)) else strata
  # 
  # # if(length(strata)==0){
  # #Calculated unique time and ties for the baseline calculation:
  # unique_time   <- unique(time)
  # 
  # tieseq <- NULL
  # index  <- NULL
  # for (i in 1:length(unique(time))) {
  #   tieseq[i] <- length(which(time==unique(time)[i]))
  #   index[[i]]  <- (which(time==unique(time)[i]))
  # }
  # 
  # #call Rcpp
  # theta_IC <- fit$ctrl.pts
  # B.spline=splines::bs(unique_time,knots=fit$internal.knots,intercept=TRUE,degree=3)
  # k=ncol(fit$bases)
  # result1 <- Lambda_estimate_ties2(knot = k, delta = event,
  #                                  z = as.matrix(z), b_spline = as.matrix(B.spline),
  #                                  theta = theta_IC, tieseq = tieseq)
  # lambda   <- result1$lambda
  # time2 <- unique(time)
  # Lambda <- cumsum(lambda)
  # 
  # baselinedata <- data.frame(unique(time),lambda,Lambda)
  # # colnames(baselinedata) <- c("time", "hazard", "Lambda")
  # 
  # res <- list("time" = unique(time),
  #             "hazard" = as.numeric(lambda),
  #             "cumulHaz" = Lambda)
  # 
  # class(res) <- "baseline"
  # 
  # return(res)
  
}






#' plotting the baseline hazard
#' 
#' Plotting the baseline hazard from a fitted `baseline` object.
#'
#' @param fit fitted object from `baseline` function.
#' @param xlab the title for the x axis.
#' @param ylab the title for the y axis.
#' @param xlim the limits of the x axis.
#' @param ylim the limits of the y axis.
#' @param title the title for the plot.

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
#' plot(base.est)
#' }
plot.baseline <- function(fit, xlab, ylab, xlim, ylim, title){
  
  
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
