#' test the proportional hazards assumption from a `coxtv` object
#' 
#' @param fit fitted \code{"coxtv"} or \code{"coxtp"}  model
#' @param parm the names of parameter to be tested
#' 
#' 
#' @examples 
#' data(ExampleData)
#' z <- ExampleData$x
#' time <- ExampleData$time
#' event <- ExampleData$event
#' fit <- coxtv(event = event, z = z, time = time)
#' tvef.ph(fit)
#' 
#' @export
tvef.ph <- function(fit, parm) {
  if (missing(fit)) stop ("Argument fit is required!")
  if (class(fit)!="coxtv" & class(fit)!="coxtp") stop("Object fit is not of class 'coxtv' or 'coxtp!")
  nsplines <- attr(fit, "nsplines"); spline <- attr(fit, "spline")
  term.ti <- names(fit$tief); term.tv <- rownames(fit$ctrl.pts)
  method <- attr(fit,"control")$method
  
  if (missing(parm)) {
    parm <- term.tv
  } else if (length(parm)>0) {
    indx <- pmatch(parm, term.tv, nomatch=0L)
    if (any(indx==0L))
      stop(gettextf("%s not matched!", parm[indx==0L]), domain=NA)
  } else stop("Invalid parm!")
  
  rownames.info <- c(rep(term.tv, each=nsplines), term.ti)
  
  if (method=="Newton") {
    invinfo <- solve(fit$info)
  } else if (method=="ProxN") {
    invinfo <- solve(fit$info+diag(sqrt(.Machine$double.eps),dim(fit$info)[1]))
  }
  # if (spline=="B-spline") {
    mat.contrast <- diff(diag(nsplines))
    ctrl.pts <- matrix(fit$ctrl.pts[term.tv%in%parm,], ncol=nsplines)
    mat.test <- sapply(parm, function(tv) {
      bread <- mat.contrast%*%ctrl.pts[parm%in%tv,]
      idx <- rownames.info%in%tv
      meat <- solve(mat.contrast%*%invinfo[idx,idx]%*%t(mat.contrast))
      stat <- t(bread)%*%meat%*%bread
      p.value <- pchisq(stat, nsplines-1, lower.tail=F)
      return(c(stat, nsplines-1, p.value))})
    colnames(mat.test) <- parm
    rownames(mat.test) <- c("chisq", "df", "p")
    return(t(mat.test))
  # } else if (spline=="P-spline") {
    
  # }
}




#' test the significance of the covariates from a `coxtv` or `coxtp` object
#' 
#' @param fit fitted \code{"coxtv"} or \code{"coxtp"}  model
#' @param parm the names of parameter to be tested
#' @param times the time points to test if the covariate is significant or not
#' 
#' @examples 
#' data(ExampleData)
#' z <- ExampleData$x
#' time <- ExampleData$time
#' event <- ExampleData$event
#' fit <- coxtv(event = event, z = z, time = time)
#' tvef.zero(fit)
#' 
#' @export
tvef.zero <- function(fit, times, parm) {
  if (missing(fit)) stop ("Argument fit is required!")
  if (class(fit)!="coxtv" & class(fit)!="coxtp") stop("Object fit is not of class 'coxtv' or 'coxtp!")
  if (missing(times)) {
    times <- fit$times
  } else {
    if (!is.numeric(times) | min(times)<0) stop("Invalid times!")
  }
  term.ti <- names(fit$tief); term.tv <- rownames(fit$ctrl.pts)
  spline <- attr(fit, "spline"); nsplines <- attr(fit, "nsplines")
  degree <- attr(fit, "degree"); knots <- fit$internal.knots
  method <- attr(fit,"control")$method
  if (missing(parm)) {
    parm <- term.tv
  } else if (length(parm)>0) {
    indx <- pmatch(parm, term.tv, nomatch=0L)
    if (any(indx==0L))
      stop(gettextf("%s not matched!", parm[indx==0L]), domain=NA)
  }
  rownames.info <- c(rep(term.tv, each=nsplines), term.ti)
  if (method=="Newton") {
    invinfo <- solve(fit$info)
  } else if (method=="ProxN") {
    invinfo <- solve(fit$info+diag(sqrt(.Machine$double.eps),dim(fit$info)[1]))
  }
  # if (spline=="B-spline") {
    bases <- splines::bs(times, degree=degree, intercept=T, knots=knots, 
                         Boundary.knots=range(fit$times))
    ctrl.pts <- matrix(fit$ctrl.pts[term.tv%in%parm,], ncol=nsplines)
    ls <- lapply(parm, function(tv) {
      est <- bases%*%ctrl.pts[parm%in%tv,]
      se <- sqrt(apply(bases, 1, function(r) {
        idx <- rownames.info%in%tv
        return(t(r)%*%invinfo[idx, idx]%*%r)}))
      stat <- est / se
      p.upper <- pnorm(stat, lower.tail=F)
      p.value <- 2*pmin(p.upper, 1-p.upper)
      mat <- cbind(est, se, stat, p.value)
      colnames(mat) <- c("est", "se", "z", "p")
      rownames(mat) <- times
      return(mat)})
    names(ls) <- parm
  # } else if (spline=="P-spline") {
  #   
  # }
  return(ls)
}

