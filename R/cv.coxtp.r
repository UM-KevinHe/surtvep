#' fit a cox Non-proportional Hazards model with P-spline or Smoothing-spline, penalization coefficient is provided by cross validation
#' 
#' Fit a cox Non-proportional Hazards model via penalized maximum likelihood. 
#' 
#'
#' @param event response variable of length `nobs`, should be a vector containing 0 or 1
#' @param z input covariate matrix, of dimension `nobs` * `nvars`; each row is an observation vector. 
#' @param time follow up time, should be a vector with non-negative numeric value
#' @param strata stratification group defined in the data. If there exists stratification group, please enter as vector. 
#' By default a non-stratified model would be implemented
#' @param penalty a character string specifying the spline term for Penalized Newton's Method. 
#' This term is added to the log-partial likelihood as the new objective function to control the model's smoothness. 
#' 
#' `P-spline` stands for "Penalized B-spline". It combines the B-spline basis with a discrete quadratic penalty on the basis coefficients between adjacent knots.
#' 
#' `Smooth-spline` refers to the Smoothing-spline, the derivative based penalties combined with B-splines. Default value is `Smooth-spline`.
#' 
#' @param lambda a user specified `lambda.spline` sequence as the penalization coefficients in front of the spline term specified by `spline`. 
#' This is the tuning parameter for penalization. Users can use `IC` to select the best tuning parameter based on the information criteria. 
#' Users can specify for larger values when the estimated time-varying effects are too high.
#' Default is `0` which refers to Newton's Method without penalization. 
#' 
#' @param nsplines number of basis functions in the B-splines to span the time-varying effects, default value is 8. 
#' We use the r function `splines::bs` to generate the B-splines. 
#' @param knots the internal knot locations (breakpoints) that define the B-splines.
#' The number of the internal knots should be \eqn{`nsplines`-`degree`-1}.
#' If `NULL`, the locations of knots are chosen to include an equal number of events within each time interval. This leads to more stable results in most cases.
#' Users can specify the internal knot locations by themselves.
#' @param degree degree of the piecewise polynomial for generating the B-spline basis functions---default is 3 for cubic splines. 
#' `degree = 2` results in the quadratic B-spline basis functions.
#' 
#' If `penalty` is `Smooth-spline`, different choices of `degree` give different results.
#' When `degree=3`, we use the cubic B-spline penalizing the second order derivative, which reduces to a linear term.
#' When `degree=2`, we use the quadratic B-spline penalizing first order derivative, which reduces to a constant. See Wood (2016) for details.
#' Default is `degree=2`.
#' 
#' @param ties a character string specifying the method for tie handling. If there are no tied
#' death times, the methods are equivalent.  By default **`ties="Breslow"`** uses the Breslow approximatio, this can be faster when many ties occured.
#' @param stop a character string specifying the stopping rule to determine convergence. Use `loglik(m)` to denote the log-partial likelihood at iteration step m.  
#' `"incre"` means we stop the algorithm when Newton's increment is less than the `tol`.
#' `"relch"` means we stop the algorithm when the `loglik(m)` divided by the  `loglik(0)` is less than the `tol`.
#' `"ratch"` means we stop the algorithm when `(loglik(m)-loglik(m-1))/(loglik(m)-loglik(0))` is less than the `tol`.
#' `"all"` means we stop the algorithm when all the stopping rules `"incre"`, `"relch"` and `"ratch"` is met. Default value is `ratch`.
#' @param tol convergence threshold for Newton's method. The algorithm continues until the method selected using `stop` converges.
#'  The default value is  `1e-6`.
#' @param iter.max maximum Iteration number, default value is  `20L`.
#' @param method a character string specifying whether to use Newton's method or Proximal Newton's method.  If `"Newton"` then exact hessian is used, 
#' while default method `"ProxN"` implementing the proximal method which can be faster and more stable when there exists ill-conditioned second-order information of the log-partial likelihood.
#' See details in Wu et al. (2022).
#' 
#' @param gamma parameter for Proximal Newton's Method `"ProxN"`. Default value is `1e8`.
#' @param btr a character string specifying the backtracking line search approach. `"dynamic"` is typical way to perform backtracking linesearch. 
#' `"static"` limits the Newton's increment, and can achieve more stable results in some extreme cases such as ill-conditioned second-order information of the log-partial likelihood, 
#' which usually occur when some predictors are categorical with low frequency for some categories. 
#' Users should be careful with `static` as this may lead to underfitting.
#' @param tau a scalar in (0,1) used to control the step size inside the back tracking linesearch. Default value is `0.5`.
#' @param parallel if `TRUE` then parallel computation is enabled. The number of threads to be used is determined by `threads`.
#' @param threads an integer indicating the number of threads to be used for parallel computation. Default is `2`. If `parallel` is false, then the value of `threads` has no effect.
#' @param fixedstep if `TRUE`, the algorithm will be forced to run `iter.max` steps regardless of the stopping criterion specified.
#' @param nfolds number of folds - default is 5. Although nfolds can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable is nfolds=3
#' @param foldid an optional vector of values between 1 and nfold identifying what fold each observation is in. If supplied, nfold can be missing.
#'
#' @return
#' \item{model.cv}{a \code{"coxtp"} object with penalization coefficient chosen based on cross validation} 
#' \item{lambda.spline}{the values of `lambda.spline` used in the fits.}
#' \item{cv.loglik}{The mean log-partial likelihood on each fold after based on the model trained on remaining folds - a vector of length length(lambda.spline).}
#' \item{lambda.selected}{Value of `lambda.spline` that gives maximum cv.loglik}
#' @export
#' 
#' @details The function runs `coxtp` `nfolds` times; each is to compute the fit with each of the folds omitted.
#'
#' @examples 
#' data(ExampleData)
#' z <- ExampleData$x
#' time  <- ExampleData$time
#' event <- ExampleData$event
#' 
#' lambda.spline  = c(0,1)
#' fit   <- cv.coxtp(event = event, z = z, time = time, lambda=lambda.spline, nfolds = 5)
cv.coxtp <- function(event , z, time, strata=NULL, nfolds = 5, foldid = NULL, penalty="Smooth-spline", nsplines=8, ties="Breslow",
                     tol=1e-9, iter.max=20L, method="ProxN", lambda=1e8,
                     btr="dynamic", tau=0.5,
                     stop="ratch", parallel=FALSE, threads=1L, degree=3L, 
                     lambda.spline = 0, fixedstep = FALSE){
  
  lambda.spline = lambda
  spline = penalty
  ord = degree + 1
  TIC_prox = FALSE
  ICLastOnly = FALSE
  
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
