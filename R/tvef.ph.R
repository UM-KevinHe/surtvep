#' get confidence interval from a `coxtv` object
#' 
#' @param fit fitted \code{"coxtv"} model
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
  if (class(fit)!="coxtv") stop("Object fit is not of class 'coxtv'!")
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
    invinfo <- solve(fit$VarianceMatrix)
  } else if (method=="ProxN") {
    invinfo <- solve(fit$VarianceMatrix+diag(sqrt(.Machine$double.eps),dim(fit$VarianceMatrix)[1]))
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
