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
#' \item{time}{the unique observed failure times.} 
#' \item{hazard}{the baseline hazard corresponding to each unique failure time point.}
#' \item{cumulHaz}{the cumulative baseline hazard corresponding to each unique failure time point.}
#'
#' @examples
#' data(ExampleData)
#' z <- ExampleData$z
#' time  <- ExampleData$time
#' event <- ExampleData$event
#' 
#' fit   <- coxtv(event = event, z = z, time = time)
#' base.est <- baseline(fit)
#' 
#' 
baseline <-  function(fit){

  if (missing(fit)) stop ("Argument fit is required!")
  if (!inherits(fit,"coxtp") & !inherits(fit,"coxtv")) stop("Object fit is not of class 'coxtv' or 'coxtp'!")
  
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
}






#' plotting the baseline hazard
#' 
#' Plotting the baseline hazard from a fitted `baseline` object.
#'
#' @param x fitted object from `baseline` function.
#' @param xlab the title for the x axis.
#' @param ylab the title for the y axis.
#' @param xlim the limits of the x axis.
#' @param ylim the limits of the y axis.
#' @param title the title for the plot.
#' @param \dots Other graphical parameters to plot
#' 
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon theme_bw theme element_text element_blank margin labs ggtitle theme_classic
#' @importFrom tibble tibble
#' @importFrom rlang .data
#' @exportS3Method plot baseline
#'
#' @examples
#' data(ExampleData)
#' z <- ExampleData$z
#' time  <- ExampleData$time
#' event <- ExampleData$event
#' 
#' fit   <- coxtv(event = event, z = z, time = time)
#' base.est <- baseline(fit)
#' plot(base.est)
plot.baseline <- function(x, xlab, ylab, xlim, ylim, title, ...){
  
  
  if (missing(xlab)) xlab <- "time"
  if (missing(ylab)) ylab <- "cumulative hazard"
  
  if (missing(x)) stop ("Argument x is required!")
  fit <- x
  if (!inherits(fit,"baseline")) stop("Object fit is not of class 'baseline'!")
  
  missingxlim <- missing(xlim); missingylim <- missing(ylim); 
  missingtitle <- missing(title);
  
  plot_data <- tibble(
    time = fit$time,
    cumulHaz = fit$cumulHaz
  )
  
  plt <- ggplot(data = plot_data, aes(x = .data$time, y = .data$cumulHaz)) +
    geom_line(size = 0.9)+
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
