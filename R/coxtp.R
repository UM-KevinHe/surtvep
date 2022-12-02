#' fit a cox Non-proportional Hazards model with P-spline or Smoothing-spline. Penalization tuning parameter
#' 
#' Fit a cox Non-proportional Hazards model via penalized maximum likelihood. 
#' 
#'
#' @param event failure events response variable of length `nobs`, where `nobs` denotes the number of observations. should be a vector containing 0 or 1
#' @param z input covariate matrix, of dimension `nobs` \code{x} `nvars`; each row is an observation vector. 
#' @param time obeserved event time, should be a vector with non-negative numeric value
#' @param strata stratification group defined in the data used for stratified model. If there exists stratification group, please enter as vector. 
#' By default a non-stratified model would be implemented
#' @param penalty a character string specifying the spline term for Penalized Newton's Method. 
#' This term is added to the log-partial likelihood as the new objective function to control the smoothness of the time-varying covariates.
#' Default is `P-spline`. Three options are `P-spline`, `Smooth-spline` and `NULL`. If `NULL`, the method will be the same as `coxtv` and `lambda` 
#' will be set as 0.
#' 
#' `P-spline` stands for Penalized B-spline. It combines the B-spline basis with a discrete quadratic penalty on the difference of basis coefficients between adjacent knots. 
#' When `lambda` goes to infinity, the time-varying effects are encouraged to be constant. 
#' 
#' `Smooth-spline` refers to the Smoothing-spline, the derivative based penalties combined with B-splines. See `degree` for different choices.
#' When `degree=3`, we use the cubic B-spline penalizing the second order derivative, which reduces to a linear term when `lambda` goes to infinity.
#' When `degree=2`, we use the quadratic B-spline penalizing first order derivative, which reduces to a constant when `lambda` goes to infinity. See Wood (2016) for details.
#' 
#' If `P-spline` or `Smooth-spline`, then `lambda` is initialized as (0.1, 1, 10). Users can modify `lambda`. See details in `lambda`.
#' 
#' @param lambda a user specified `lambda.spline` sequence as the penalization coefficients in front of the spline term specified by `spline`. 
#' This is the tuning parameter for penalization. Users can use `IC` to select the best tuning parameter based on the information criteria. 
#' Users can specify for larger values when the estimated time-varying effects are too high.
#' Default is `0` which refers to Newton's Method without penalization. 
#' 
#' @param nsplines number of basis functions in the B-splines to span the time-varying effects, default value is 8. 
#' We use the r function `splines::bs` to generate the B-splines. 
#' @param knots the internal knot locations (breakpoints) that define the B-splines.
#' The number of the internal knots should be `nsplines`-`degree`-1.
#' If `NULL`, the locations of knots are chosen to include an equal number of events within each time interval. This leads to more stable results in most cases.
#' Users can specify the internal knot locations by themselves.
#' @param degree degree of the piecewise polynomial for generating the B-spline basis functions---default is 3 for cubic splines. 
#' `degree = 2` results in the quadratic B-spline basis functions.
#' 
#' If `penalty` is `Smooth-spline`, different choices of `degree` give different results.
#' When `degree=3`, we use the cubic B-spline penalizing the second order derivative, which reduces to a linear term when `lambda` goes to infinity.
#' When `degree=2`, we use the quadratic B-spline penalizing first order derivative, which reduces to a constant when `lambda` goes to infinity. See Wood (2016) for details.
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
#'
#' @return An object with S3 class \code{"coxtp"}. 
#' \item{call}{the call that produced this object}
#' \item{beta}{estimated coefficient matrix of dimension `len_unique_t` * `nvars`, where `len_unique_t` is the length of unique follow-up `time`.
#' Each row represents the coefficients at the corresponding input follow-up time}
#' \item{bases}{the basis matrix used in model fitting. If `ties="None"`, the dimension is `nvars` * `nsplines`; 
#' if `ties="Breslow"`, the dimension is `len_unique_t` * `nsplines`. The matrix is constructed using `bs::splines` function.}
#' \item{ctrl.pts}{estimated coefficient matrix of dimension `nvars` * `nsplines`. 
#' Each row represents a covariate's coefficient on the `nsplines` dimensional basis functions.} 
#' \item{Hessian}{the Hessian matrix of the log-partial likelihood, of which the dimension is `nsplines * nvars` multiplied by `nsplines * nvars`.}
#' \item{internal.knots}{the internal knot locations of the basis functions. The locations of knots are chosen to include an equal number of events within each time interval.}
#' \item{nobs}{number of observations}
#' \item{spline}{the spline type user specified.}
#' \item{theta.list}{a list of `ctrl.pts` of length `m`, contains the updated `ctrl.pts` after each algorithm iteration.}
#' \item{VarianceMatrix}{the variance matrix of the estimated function, which is the inverse of the negative `Hessian` matrix.}
#'
#'
#' @details 
#' The sequence of models implied by `lambda.spline` is fit by Newton's method (Proximal Newton's method). The objective function is
#' \deqn{loglik - P_{\lambda}},
#' where P_{lambda} can be `P-spline` or `Smooth-spline`. The \eqn{`\lambda`} is the tuning  parameter \eqn{`\lambda`}. The initial sequence can be defined by user.
#' `IC`provides different information criteria to choose the tuning parameter \eqn{`\lambda`}. `cv.coxtp` uses  the cross validation to choose the tuning parameter.
#'
#' @seealso \code{coef}, \code{plot}, \code{IC} and \code{cv.coxtp}.
#' 
#' 
#' @export
#'
#' @examples 
#' data(ExampleData)
#' z <- ExampleData$x
#' time  <- ExampleData$time
#' event <- ExampleData$event
#' 
#' lambda  = c(0,1)
#' fit   <- coxtp(event = event, z = z, time = time, lambda=lambda)
#' 
#' 
#' 
#' @references 
#' Gray, R.~J.
#' \emph{Flexible methods for analyzing survival data using splines, with applications to breast cancer prognosis. (1992), Journal of the American Statistical Association, Vol. 87, 942--951}.
#' \cr
#' 
#' Gray, R.~J.
#' \emph{Spline-based tests in survival analysis. (1994), Biometrics, Vol. 50, 640--652}.
#' \cr
#' 
#' Lingfeng Luo, Kevin He, Wenbo Wu and Jeremy M.G. Taylor 
#' \emph{Using Information Criteria to Select Smoothing Parameters when Analyzing Survival Data with Time-Varying Coefficient Hazard Models (2022)}.
#' \cr
#' 
#' Wenbo Wu, Jeremy M G Taylor, Andrew F Brouwer, Lingfeng Luo, Jian Kang, Hui Jiang and Kevin He. 
#' \emph{Scalable proximal Methods for cause-specific hazard modeling with time-varying coefficients (2022), Lifetime Data Analysis, Vol. 28(2), 194-218}.
#' \cr
#' 
#' Wood, Simon N.
#' \emph{P-splines with derivative based penalties and tensor product smoothing of unevenly distributed data. (2017) 
#' Statistics and Computing, Vol. 27(4), 985-989}.
#' 
#' 
#' 
#' 
coxtp <- function(event , z , time ,strata=NULL ,penalty="Smooth-spline", nsplines=8, 
                  lambda = 0, degree=3L,
                  knots = NULL,
                  ties="Breslow",
                  tol=1e-9, iter.max=20L, method="ProxN", gamma=1e8,
                  btr="dynamic", tau=0.5,
                  stop="ratch", parallel=FALSE, threads=1L, 
                  fixedstep = FALSE,
                  InfoCrit = FALSE,...){
  
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
  
  InfoCrit = FALSE
  
  data_NR <- data.frame(event=event, time=time, z, strata=stratum, stringsAsFactors=F)
  Z.char <- colnames(data_NR)[c(3:(ncol(data_NR)-1))]
  fmla <- formula(paste0("Surv(time, event)~",
                         paste(c(paste0("tv(", Z.char, ")"), "strata(strata)"), collapse="+")))

  fit <- coxtp.base(fmla, data_NR, nsplines=nsplines, spline=spline, ties=ties, stop=stop,
                    method = method, btr = btr,
                    lambda_spline = lambda.spline,TIC_prox = TIC_prox, ord = ord, degree = degree,
                    tol = tol, iter.max = iter.max, tau= tau, parallel = parallel, threads = threads,
                    fixedstep = fixedstep,
                    ICLastOnly = ICLastOnly)


  return(fit)


}






#' #' fit a cox Non-proportional Hazards model with P-spline or Smoothing-spline, penalization coefficient is provided by the information criteria
#' #' 
#' #' Fit a cox Non-proportional Hazards model via penalized maximum likelihood. 
#' #' 
#' #' @details put the link to the wesite()...
#' #'
#' #' @param event response variable of length `nobs`, should be a vector containing 0 or 1
#' #' @param z input covariate matrix, of dimension `nobs` * `nvars`; each row is an observation vector. 
#' #' @param time follow up time, should be a vector with non-negative numeric value
#' #' @param strata stratification group defined in the data. If there exist stratification group, please enter as vector. 
#' #' Otherwise a non-stratified model would be implemented
#' #' @param spline a character string specifying the spline term for Penalized Newton's Method. 
#' #' This term is added to the log-partial likelihood as the new objective function to control the model's smoothness. 
#' #' `P-spline` uses P-splines, which are low rank smoothers using a B-spline basis defined on evenly spaced knots. 
#' #' `Smooth-spline` refers to the Smoothing-spline, the reduced ran spline smoothers with derivative based penalties. 
#' #' Default value is `Smooth-spline`.
#' #' @param lambda.spline A user specified `lambda` sequence as the penalization coefficients. Default is **`lambda.spline = 0`** which refers to Newton's Method without penalization.
#' #' @param nsplines Number of basis functions in the B-splines, default value is 8.  
#' #' @param ties a character string specifying the method for tie handling. If there are no tied
#' #' death times, the methods are equivalent.  **`ties="Breslow"`** uses the Breslow approximation.
#' #' @param stop a character string specifying the stopping rule to determine convergence. Use `loglik(m)` to denote the log-partial likelihood at iteration step m.  
#' #' `"incre"` means we stop the algorithm when Newton's increment is less than the `tol`.
#' #' `"relch"` means we stop the algorithm when the `loglik(m)` divided by the  `loglik(0)` is less than the `tol`.
#' #' `"ratch"` means we stop the algorithm when `(loglik(m)-loglik(m-1))/(loglik(m)-loglik(0))` is less than the `tol`.
#' #' `"all"` means we stop the algorithm when all the stopping rules `"incre"`, `"relch"` and `"ratch"` is met. Default value is `ratch`.
#' #' @param tol Convergence threshold for Newton's method. The algorithm continues until the method selected using `stop` converges.
#' #'  The default value is  `1e-6`.
#' #' @param iter.max Maximum Iteration number, default value is  `20L`
#' #' @param method a character string specifying whether to use Newton's method or Proximal Newton's method.  If `"Newton"` then exact hessian is used (default), 
#' #' while `"ProxN"` adds a small amount `lambda` to the diagonal matrix, and can be faster and more stable.
#' #' @param lambda Parameter for Proximal Newton's Method. Default is **`lambda=1e8`**
#' #' @param btr a character string specifying the backtracking line search approach. `"dynamic"` is typical way to perform backtracking linesearch. 
#' #' `"static"` limits the Newton's increment, and can achieve more stable results in some extreme cases.
#' #' @param tau a scalar in (0,1) used to control the step size inside the back tracking linesearch. Default value is `0.5`.
#' #' @param parallel If `TRUE` then parallel computation is enabled. The number of threads to be used is determined by `threads`.
#' #' @param threads an integer indicating the number of threads to be used for parallel computation. Default is `2`. If `parallel` is false, then the value of `threads` has no effect.
#' #' @param degree degree of the piecewise polynomial---default is 3 for cubic splines.
#' #' @param ord a positive integer giving the order of the spline function. 
#' #' This is the number of coefficients in each piecewise polynomial segment, thus a cubic spline has order 4. Defaults to 4.
#' #' @param fixedstep If `TRUE`, the algorithm will be forced to run `iter.max` steps regardless of the stopping criterion specified.
#' #' @param TIC_prox When calculating information criteria, there might be numerical issue(second order derivative, 
#' #' Hessian matrix approximate is singular), thus we proposed to add a small term to the diagonal. Default **`TIC_prox = FALSE`**
#' #' @param ICLastOnly Only calculate the last information criteria if is TRUE. Default is **`ICLastOnly=FALSE`**
#' #'
#' #' @return An object with S3 class \code{"coxtp"}. 
#' #' \item{call}{the call that produced this object}
#' #' \item{beta}{estimated coefficient matrix of dimension `len_unique_t` * `nvars`, where `len_unique_t` is the length of unique follow-up `time`.
#' #' Each row represents the coefficients at the corresponding input follow-up time}
#' #' \item{bases}{the basis matrix used in model fitting. If `ties="None"`, the dimension is `nvars` * `nsplines`; 
#' #' if `ties="Breslow"`, the dimension is `len_unique_t` * `nsplines`. The matrix is constructed using `bs::splines` function.}
#' #' \item{ctrl.pts}{estimated coefficient matrix of dimension `nvars` * `nsplines`. 
#' #' Each row represents a covariate's coefficient on the `nsplines` dimensional basis functions.} 
#' #' \item{Hessian}{the Hessian matrix of the log-partial likelihood, of which the dimension is `nsplines * nvars` multiplied by `nsplines * nvars`.}
#' #' \item{internal.knots}{the internal knot locations of the basis functions. The locations of knots are chosen to include an equal number of events within each time interval.}
#' #' \item{nobs}{number of observations}
#' #' \item{spline}{The spline type user specified.}
#' #' \item{theta.list}{a list of `ctrl.pts` of length `m`, contains the updated `ctrl.pts` after each algorithm iteration.}
#' #' \item{VarianceMatrix}{The variance matrix of the estimated function, which is the inverse of the negative `Hessian` matrix.}
#' #'
#' #' @export
#' #'
#' #' @examples 
#' #' data(ExampleData)
#' #' z <- ExampleData$x
#' #' time  <- ExampleData$time
#' #' event <- ExampleData$event
#' #' 
#' #' lambda.spline  = c(0,1)
#' #' fit   <- ic.coxtp(event = event, z = z, time = time, lambda.spline=lambda.spline)
#' ic.coxtp <- function(event , z , time ,strata=c() ,spline="Smooth-spline", nsplines=8, ties="Breslow",
#'                   tol=1e-9, iter.max=20L, method="ProxN", lambda=1e8,
#'                   btr="dynamic", tau=0.5,
#'                   stop="ratch", parallel=FALSE, threads=1L, degree=3L, TIC_prox = FALSE,
#'                   lambda.spline = 0, ord = 4, fixedstep = FALSE){
#'   
#'   if(length(strata)==0){
#'     stratum=rep(1, length(time))
#'   } else {
#'     stratum=strata
#'   }
#'   
#'   InfoCrit = TRUE
#'   
#'   data_NR <- data.frame(event=event, time=time, z, strata=stratum, stringsAsFactors=F)
#'   Z.char <- colnames(data_NR)[c(3:(ncol(data_NR)-1))]
#'   fmla <- formula(paste0("Surv(time, event)~",
#'                          paste(c(paste0("tv(", Z.char, ")"), "strata(strata)"), collapse="+")))
#'   
#'   lambda_all <- lambda.spline
#'   model <- list()
#'   
#'   for(lambda_index in 1:length(lambda_all)){
#'     
#'     model1 <- coxtp.base(fmla, data_NR, nsplines=nsplines, spline=spline, ties=ties, stop=stop,
#'                          method = method, btr = btr,
#'                          lambda_spline = lambda_all[lambda_index],TIC_prox = TIC_prox, ord = ord, degree = degree,
#'                          tol = tol, iter.max = iter.max, tau= tau, parallel = parallel, threads = threads,
#'                          fixedstep = fixedstep,
#'                          ICLastOnly = InfoCrit)
#'     
#'     model[[lambda_index]] <- model1
#'   }
#'   
#'   
#'   AIC_all <- NULL
#'   TIC_all <- NULL
#'   GIC_all <- NULL
#'   
#'   for (i in 1:length(lambda_all)){
#'     AIC_all<-c(AIC_all, model[[i]]$AIC_all)
#'     TIC_all<-c(TIC_all, model[[i]]$TIC_all)
#'     GIC_all<-c(GIC_all, model[[i]]$GIC_all)
#'   }
#'   
#'   AIC_lambda <- which.min(AIC_all)
#'   TIC_lambda <- which.min(TIC_all)
#'   GIC_lambda <- which.min(GIC_all)
#'   
#'   model_AIC <- model[[AIC_lambda]]
#'   model_TIC <- model[[TIC_lambda]]
#'   model_GIC <- model[[GIC_lambda]]
#'   
#'   lambda.selected <- data.frame(c(lambda_all[AIC_lambda], lambda_all[TIC_lambda], lambda_all[GIC_lambda]))
#'   rownames(lambda.selected) <- c("AIC","TIC","GIC")
#'   colnames(lambda.selected) <- c("value")
#'   
#'   z_names=colnames(z)
#'   p        <- ncol(z)
#'   
#'   
#'   res = NULL
#'   res$model.AIC = model_AIC
#'   res$model.TIC = model_TIC
#'   res$model.GIC = model_GIC
#'   res$lambda.selected = lambda.selected
#'   # res$p=p
#'   # res$z_names=z_names
#'   
#'   class(res) <- "ic.coxtp"
#'   
#'   return(res)
#'   
#'   
#' }
