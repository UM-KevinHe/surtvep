#' helper function to get time-varying coefficients
#' 
#' Gives the time-varying coefficients based on a fitted `coxtv` or `coxtp` subject. Users can specify  the observation  time.
#'
#' @param fit model from `coxtv` or `coxtp`.
#' @param times the observation time. If `NULL`, the observation time for fitting the model will be used.
#' 
#' @return a matrix of the time-varying coefficients. The dimension is `length(times)` x `nvars`, where `nvars` is the number
#' of covariates in the fitted mode.
#' Each row represents the time-varying coefficients at the corresponding time.
#' ``
#' @export
#' 
#' @examples 
#' z     <- ExampleData$x
#' time  <- ExampleData$time
#' event <- ExampleData$event
#' fit   <- coxtv(event = event, z = z, time = time, degree = 2)
#' coef  <- get.tvcoef(fit)
#' 
#' 
get.tvcoef <- function(fit, times) {
  if (missing(fit)) stop ("Argument fit is required!")
  if (!class(fit)%in%c("coxtv","coxtp")) stop("Object fit is not of the classes 'coxtv' or 'coxtp'!")
  if (missing(times)) times <- fit$times
  if (!is.numeric(times) | min(times)<0) stop("Invalid times!")
  times <- times[order(times)]; nsplines <- attr(fit, "nsplines")
  spline <- attr(fit, "spline"); degree <- attr(fit, "degree.spline")
  knots <- fit$internal.knots; term.tv <- rownames(fit$ctrl.pts)
  # if (missing(parm)) {
  parm <- term.tv
  # } else if (length(parm)>0) {
  # indx <- pmatch(parm, term.tv, nomatch=0L)
  # if (any(indx==0L))
  # stop(gettextf("%s not matched!", parm[indx==0L]), domain=NA)
  # } else stop("Invalid parm!")
  # if (spline=="B-spline") {
  bases <- splines::bs(times, degree=degree, intercept=T, knots=knots, 
                       Boundary.knots=range(fit$times))
  # int.bases <- splines2::ibs(times, degree=degree, intercept=T, knots=knots, 
  #                            Boundary.knots=range(fit$times))
  ctrl.pts <- matrix(fit$ctrl.pts[term.tv%in%parm,], ncol=nsplines)
  mat.tvef <- bases%*%t(ctrl.pts) 
  # mat.cumtvef <- int.bases%*%t(ctrl.pts)
  colnames(mat.tvef) <- parm 
  # colnames(mat.cumtvef) <- parm
  rownames(mat.tvef) <- times
  # rownames(mat.cumtvef) <- times
  # ls <- list(tvef=mat.tvef)
  return(mat.tvef)
  # return(ls)
  # } else if (spline=="P-spline") {
  
  # }
}



