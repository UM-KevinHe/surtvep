#' fit a cox Non-proportional Hazards model with P-spline or Smoothing-spline, penalization coefficient is provided by cross validation
#' 
#' Fit a cox Non-proportional Hazards model via penalized maximum likelihood. 
#' 
#'
#' @param event event vector, should be a vector containing 0 or 1
#' @param z Covariate matrix
#' @param time Time vector, should be a vector with non-negative numeric value
#' @param strata stratification group defined in the data. If there exist stratification group, please enter as vector.
#' @param spline The spline term for Penalized Newton's Method(Add section 
#' Number related to the paper). Default setting is **`spline="Smooth-spline"`**
#' @param nsplines Number of base functions in the B-splines, default is 8.
#' @param ties Ways to deal with ties, default is **`ties="Breslow"`**:
#' @param tol Convergence threshold. The default threshold is set as **`tol=1e-6`**
#' @param iter.max Maximum Iteration number, default is **`iter.max=20L`**
#' @param method Selecting Method used, default is **`method="ProxN"`**
#' @param lambda Parameter for Proximal Newton's Method. Default is **`lambda=1e8`**
#' @param btr Backtracking line search approach, default is **`btr="dynamic"`**:
#' @param tau In the Newton's Method, Default is **`tau=0.5`**. Used to control for step size.
#' @param stop Stopping rule, default is **`stop="ratch"`**.
#' @param parallel Parallel computation, Default is **`parallel=FALSE`**
#' @param threads Parallel computation parameter(number of cores)Default is **`threads=1L`**
#' @param degree Degree of smoothing spline. Default setting is **`degree=3L`**.
#' @param TIC_prox When calculating information criteria, there might be numerical issue(second order derivative, 
#' Hessian matrix approximate is singular), thus we proposed to add a small term to the diagonal. Default **`TIC_prox = FALSE`**
#' @param lambda.spline  Smoothing parameter lambda. Default is **`lambda.spline = 0`** which refers to Newton's Method without penalization.
#' @param ord Specify which derivative to penalize. Default setting is **`ord=4`**.
#' @param fixedstep There might be times when the stopping criteria not working, thus, 
#' the number of steps could be set manually. Default value is **`fixedstep = FALSE`**, if it is true, will stop by `iter.max`
#' @param ICLastOnly Only calculate the last information criteria if is TRUE. Default is **`ICLastOnly=FALSE`**
#'
#' @return
#' \item{model.cv}{a \code{"coxtp"} object with penalization coefficient chosen based on cross validation} 
#' \item{lambda.selected}{The selected penalization coefficients for different information criteria}
#' @export
#'
#' @examples 
#' data(ExampleData)
#' z <- ExampleData$x
#' time  <- ExampleData$time
#' event <- ExampleData$event
#' 
#' lambda.spline  = c(0,1)
#' fit   <- ic.coxtp(event = event, z = z, time = time, lambda.spline=lambda.spline)
cv.coxtp <- function(event , z , time ,strata=c() ,spline="Smooth-spline", nsplines=8, ties="Breslow",
                     tol=1e-9, iter.max=20L, method="ProxN", lambda=1e8,
                     btr="dynamic", tau=0.5,
                     stop="ratch", parallel=FALSE, threads=1L, degree=3L, TIC_prox = FALSE,
                     lambda.spline = 0, ord = 4, fixedstep = FALSE){
  
  if(length(strata)==0){
    stratum=rep(1, length(time))
  } else {
    stratum=strata
  }
  
  InfoCrit = TRUE
  
  data_NR <- data.frame(event=event, time=time, z, strata=stratum, stringsAsFactors=F)
  Z.char <- colnames(data_NR)[c(3:(ncol(data_NR)-1))]
  fmla <- formula(paste0("Surv(time, event)~",
                         paste(c(paste0("tv(", Z.char, ")"), "strata(strata)"), collapse="+")))
  
  lambda_all <- lambda.spline
  model <- list()
  
  for(lambda_index in 1:length(lambda_all)){
    
    model1 <- coxtp.base(fmla, data_NR, nsplines=nsplines, spline=spline, ties=ties, stop=stop,
                         method = method, btr = btr,
                         lambda_spline = lambda_all[lambda_index],TIC_prox = TIC_prox, ord = ord, degree = degree,
                         tol = tol, iter.max = iter.max, tau= tau, parallel = parallel, threads = threads,
                         fixedstep = fixedstep,
                         ICLastOnly = InfoCrit)
    
    model[[lambda_index]] <- model1
  }
  
  
  res = NULL
  res$model.cv = model1
  res$lambda.selected = lambda_index

  
  class(res) <- "cv.coxtp"
  
  return(res)
  
  
}
