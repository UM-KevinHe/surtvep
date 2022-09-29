
#' Cox Non-proportional Hazards model:
#' 
#' Descritpion (...)
#' 
#' @details put the link to the wesite()...
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
#' @param method Selecting Method used, default is **`method="Newton"`**
#' @param lambda Parameter for Proximal Newton's Method. Default is **`lambda=1e8`**
#' @param btr Backtracking line search approach, default is `btr="static"`:
#' @param tau (Alpha or beta?) in the Newton's Method, Default is **`tau=0.5`**. Used to control for step size.
#' @param stop Stopping rule, default is **`stop="ratch"`**:
#' @param parallel Parallel computation, Default is **`parallel=FALSE`**
#' @param threads Parallel computation parameter(number of cores)Default is **`threads=1L`**
#' @param degree Degree of smoothing spline. Default setting is **`degree=3L`**.
#' @param ord Specify which derivative to penalize. Default setting is **`ord=4`**.
#' @param fixedstep There might be times when the stopping criteria not working, thus, 
#' the number of steps could be set manually. Default value is **`fixedstep = FALSE`**, if it is true, will stop by `iter.max`
#'
#' @return
#' @export
#'
#' @examples 
coxtp <- function(event , z , time ,strata=c() ,spline="Smooth-spline", nsplines=8, ties="Breslow",
                    tol=1e-9, iter.max=20L, method="Newton", lambda=1e8,
                    btr="static", tau=0.5,
                    stop="ratch", parallel=FALSE, threads=1L, degree=3L,
                    lambda_spline = 0, ord = 4, fixedstep = FALSE){
  
  if(length(strata)==0){
    stratum=rep(1, length(time))
  } else {
    stratum=strata
  }

  data_NR <- data.frame(event=event, time=time, z, strata=stratum, stringsAsFactors=F)
  #Z.char <- paste0("X", 1:p)
  Z.char <- colnames(data_NR)[c(3:(ncol(data_NR)-1))]
  fmla <- formula(paste0("Surv(time, event)~",
                         paste(c(paste0("tv(", Z.char, ")"), "strata(strata)"), collapse="+")))


  model1 <- coxtv.base(fmla, data_NR, nsplines=nsplines, spline=spline, ties=ties, stop=stop,
                  method = method, btr = btr,
                  lambda_spline = lambda_all[lambda_index],TIC_prox = TIC_prox, ord = ord, degree = degree,
                  tol = tol, iter.max = iter.max, tau= tau, parallel = parallel, threads = threads,
                  fixedstep = fixedstep)


}

