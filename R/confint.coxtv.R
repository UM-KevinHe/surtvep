#' get confidence intervals of time-varying coefficients from a fitted object
#' 
#' Get confidence intervals of time-varying coefficients from a fitted `coxtv` or `coxtp` object. 
#' 
#' @param object fitted \code{"coxtv"} model.
#' @param parm the names of parameters.
#' @param level the confidence level. The default value is 0.95.
#' @param time the time points for which the confidence intervals to be estimated. 
#' The default value is the unique observed event times in the dataset fitting the time-varying effects model.
#' @param \dots other parameters to function
#' 
#' @return A list where each element corresponds to one of the parameters specified in `parm`. Each element in the 
#' list is a matrix, with rows corresponding to the specified `time` points and three columns representing the 
#' estimated values of the parameter, and the lower and upper bounds of the confidence interval at the specified 
#' confidence `level`. The length of the list is determined by the number of parameters in `parm`, and each matrix 
#' has rows equal to the number of specified `time` points.
#' 
#' @examples 
#' data(ExampleData)
#' z <- ExampleData$z
#' time <- ExampleData$time
#' event <- ExampleData$event
#' fit <- coxtv(event = event, z = z, time = time)
#' confint <- confint(fit)
#' 
#' @exportS3Method confint coxtv
confint.coxtv <- function(object, parm, level=0.95, time, ...) {
  if (missing(object)) stop ("Argument object is required!")
  if (!inherits(object,"coxtv")) stop("object is not of class 'coxtv'!")
  fit <- object
  if (missing(time)) {
    time <- fit$times
  } else {
    if (!is.numeric(time) | min(time)<0) stop("Invalid time!")
  }
  if (!is.numeric(level) | level[1]>1 | level[1]<0) stop("Invalid level!")
  level <- level[1]
  time <- time[order(time)]
  time <- unique(time)
  spline <- attr(fit, "spline"); degree <- attr(fit, "degree.spline")
  knots <- fit$internal.knots; nsplines <- attr(fit, "nsplines")
  method <- attr(fit, "control")$method
  term.tv <- rownames(fit$ctrl.pts)
  if (missing(parm)) {
    parm <- term.tv
  } else if (length(parm)>0) {
    indx <- pmatch(parm, term.tv, nomatch=0L)
    if (any(indx==0L))
      stop(gettextf("%s not matched!", parm[indx==0L]), domain=NA)
  } else stop("Invalid parm!")
  rownames.info <- rep(term.tv, each=nsplines)
  if (method=="Newton") {
    invinfo <- solve(fit$VarianceMatrix)
  } else if (method=="ProxN") {
    invinfo <- solve(fit$VarianceMatrix+diag(sqrt(.Machine$double.eps),dim(fit$VarianceMatrix)[1]))
  }
  # parm.ti <- intersect(parm, c(term.ti))
  parm.tv <- intersect(parm, c(term.tv))
  quant.upper <- qnorm((1+level)/2)
  ls <- list()
  # if (length(parm.ti)!=0) {
  #   est.ti <- fit$tief[term.ti%in%parm.ti]
  #   se.ti <- c(sqrt(diag(as.matrix(invinfo[rownames.info%in%parm.ti,
  #                                          rownames.info%in%parm.ti]))))
  #   mat.ti <- cbind(est.ti, est.ti-quant.upper*se.ti, est.ti+quant.upper*se.ti)
  #   colnames(mat.ti) <- 
  #     c("est", paste0(round(100*c(1-(1+level)/2,(1+level)/2),1),"%"))
  #   rownames(mat.ti) <- parm.ti
  #   ls$tief <- mat.ti
  # }
  # if (length(parm.tv)!=0) {
  # if (spline=="B-spline") {
  bases <- splines::bs(time, degree=degree, intercept=T, knots=knots, 
                       Boundary.knots=range(fit$times))
  ctrl.pts <- matrix(fit$ctrl.pts[term.tv%in%parm.tv,], ncol=nsplines)
  ls$tvef <- lapply(parm.tv, function(tv) {
    est.tv <- bases%*%ctrl.pts[parm.tv%in%tv,]
    se.tv <- sqrt(apply(bases, 1, function(r) {
      idx <- rownames.info%in%tv
      return(t(r)%*%invinfo[idx, idx]%*%r)}))
    mat.tv <- cbind(est.tv, est.tv-quant.upper*se.tv, 
                    est.tv+quant.upper*se.tv)
    colnames(mat.tv) <- 
      c("est", paste0(round(100*c(1-(1+level)/2,(1+level)/2),1),"%"))
    rownames(mat.tv) <- time
    return(mat.tv)
  })
  names(ls$tvef) <- parm.tv
  
  return(ls)
}




