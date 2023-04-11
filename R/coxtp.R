#' fit a Cox non-proportional hazards model with P-spline or Smoothing-spline, with penalization tuning parameter chosen by information criteria or cross-validation
#' 
#' Fit a Cox non-proportional hazards model via penalized maximum likelihood. 
#' 
#'
#' @param event failure events response variable of length `nobs`, where `nobs` denotes the number of observations. It should be a vector containing 0 or 1.
#' @param z input covariate matrix, with `nobs` rows and `nvars` columns; each row is an observation vector. 
#' @param time observed event time, which should be a vector with non-negative numeric values.
#' @param strata a vector of indicators defined in the data used for stratification. Such column in the data should be entered as a vector.
#' If there exists a stratification group, please enter it as a vector. 
#' By default, an unstratified model is be implemented.
#' 
#' @param penalty a character string specifying the spline term for penalized Newton method. 
#' This term is added to the log-partial likelihood, and the penalized log-partial likelihood serves as the new objective function to control the smoothness of the time-varying effects.
#' Default is `P-spline`. Three options are `P-spline`, `Smooth-spline` and `NULL`. 
#' If `NULL`, the method will be the same as `coxtv` (unpenalized time-varying effects models) and `lambda` (defined below)
#' will be set as 0. 
#' 
#' `P-spline` stands for Penalized B-spline. It combines the B-spline basis with a discrete quadratic penalty on the difference of basis coefficients between adjacent knots. 
#' When `lambda` goes to infinity, the time-varying effects are reduced to be constant. 
#' 
#' `Smooth-spline` refers to the Smoothing-spline, the derivative-based penalties combined with B-splines. See `degree` for different choices.
#' When `degree=3`, we use the cubic B-spline penalizing the second-order derivative, which reduces the time-varying effect to a linear term when `lambda` goes to infinity.
#' When `degree=2`, we use the quadratic B-spline penalizing first-order derivative, which reduces the time-varying effect to a constant when `lambda` goes to infinity. See Wood (2016) for details.
#' 
#' If `P-spline` or `Smooth-spline`, then `lambda` is initialized as a sequence (0.1, 1, 10). Users can modify `lambda`. See details in `lambda`.
#' 
#' @param lambda a user-specified `lambda` sequence as the penalization coefficients in front of the spline term specified by `penalty`. 
#' This is the tuning parameter for penalization. Users can use `IC` to select the best tuning parameter based on the information criteria. 
#' Users can specify larger values when the absolute values of the estimated time-varying effects are too large.
#' Default is `0`, which refers to Newton method without penalization. 
#' 
#' @param nsplines number of basis functions in the splines to span the time-varying effects, whose value is 8. 
#' We use the R function `splines::bs` to generate the B-splines. 
#' 
#' @param knots the internal knot locations (breakpoints) that define the B-splines.
#' The number of the internal knots should be `nsplines`-`degree`-1.
#' If `NULL`, the locations of knots are chosen as quantiles of distinct failure time points.
#' This choice leads to more stable results in most cases.
#' Users can specify the internal knot locations by themselves.
#' 
#' @param degree degree of the piecewise polynomial for generating the B-spline basis functions---default is 3 for cubic splines. 
#' `degree = 2` results in the quadratic B-spline basis functions. 
#' 
#' If the `penalty` is `P-spline` or `NULL`, `degree`'s default value is 3. 
#' 
#' If the `penalty` is `Smooth-spline`, `degree`'s default value is 2. 
#' 
#' @param ties a character string specifying the method for tie handling. If there are no tied events, 
#' the methods are equivalent.  By default `"Breslow"` uses the Breslow approximation, which can be faster when many ties are present.
#' 
#' @param stop a character string specifying the stopping rule to determine convergence. Use \eqn{loglik(m)} to denote the log-partial likelihood at iteration step m.  
#' `"incre"` means we stop the algorithm when Newton's increment is less than the `tol`.
#' `"relch"` means we stop the algorithm when the \eqn{loglik(m)} divided by the  \eqn{loglik(0)} is less than the `tol`.
#' `"ratch"` means we stop the algorithm when \eqn{(loglik(m)-loglik(m-1))/(loglik(m)-loglik(0))} is less than the `tol`.
#' `"all"` means we stop the algorithm when all the stopping rules `"incre"`, `"relch"` and `"ratch"` are met. 
#' Default value is `ratch`. If the maximum iteration steps `iter.max` is achieved, the algorithm stops before the stopping rule is met.
#' 
#' @param tol tolerance used for stopping the algorithm. The algorithm continues until the method selected using `stop` converges.
#'  The default value is  `1e-6`.
#' @param iter.max maximum iteration number if the stopping criterion specified by `stop` is not satisfied. Default value is  20.
#' @param method a character string specifying whether to use Newton method or proximal Newton method.  If `"Newton"` then Hessian is used, 
#' while the default method `"ProxN"` implements the proximal Newton which can be faster and more stable when there exists ill-conditioned second-order information of the log-partial likelihood.
#' See details in Wu et al. (2022).
#' 
#' @param gamma parameter for proximal Newton method `"ProxN"`. The default value is `1e8`.
#' @param btr a character string specifying the backtracking line-search approach. `"dynamic"` is a typical way to perform backtracking line-search. See details in Convex Optimization by Boyd and Vandenberghe (2009).
#' `"static"` limits Newton's increment and can achieve more stable results in some extreme cases, such as ill-conditioned second-order information of the log-partial likelihood, 
#' which usually occurs when some predictors are categorical with low frequency for some categories. 
#' Users should be careful with `static` as this may lead to under-fitting.
#' @param tau a scalar in (0,1) used to control the step size inside the backtracking line-search. The default value is 0.5.
#' @param parallel if `TRUE`, then the parallel computation is enabled. The number of threads in use is determined by `threads`.
#' @param threads an integer indicating the number of threads to be used for parallel computation. Default is `2`. If `parallel` is false, then the value of `threads` has no effect.
#' @param fixedstep if `TRUE`, the algorithm will be forced to run `iter.max` steps regardless of the stopping criterion specified.
#' 
#' @return A list of objects with S3 class \code{"coxtp"}. The length is the same as that of `lambda`; each represents the model output with each value of the tuning parameter `lambda`.
#' \item{call}{the call that produced this object.}
#' \item{beta}{the estimated time varying coefficient for each predictor at each unique time. It is a matrix of dimension `len_unique_t` x `nvars`, where `len_unique_t` is the length of unique follow-up `time`.
#' Each row represents the coefficients at the corresponding input observation time.}
#' 
#' \item{bases}{the basis matrix used in model fitting. If `ties="None"`, the dimension of the basis matrix is `nvars`-by-`nsplines`; 
#' if `ties="Breslow"`, the dimension is `len_unique_t`-by-`nsplines`. The matrix is constructed using the `bs::splines` function.}
#' \item{ctrl.pts}{estimated coefficient of the basis matrix of dimension `nvars`-by-`nsplines`. 
#' Each row represents a covariate's coefficient on the `nsplines` dimensional basis functions.} 
#' \item{Hessian}{the Hessian matrix of the log-partial likelihood, of which the dimension is `nsplines * nvars` -by- `nsplines * nvars`.}
#' \item{internal.knots}{the internal knot locations of the basis functions. The locations of knots are chosen to include an equal number of events within each time interval.}
#' \item{nobs}{number of observations.}
#' \item{spline}{the spline type user specified.}
#' \item{theta.list}{a list of `ctrl.pts` of length `m`, contains the updated `ctrl.pts` after each algorithm iteration.}
#' \item{VarianceMatrix}{the variance matrix of the estimated coefficients of the basis matrix, which is the inverse of the negative `Hessian` matrix.}
#'
#'
#' @details 
#' The sequence of models implied by `lambda.spline` is fit by Newton method (proximal Newton method). The objective function is
#' \deqn{loglik - P_{\lambda},}
#' where \eqn{P_{\lambda}} can be `P-spline` or `Smooth-spline`. The \eqn{\lambda} is the tuning  parameter \eqn{\lambda}. Users can define the initial sequence.
#' `IC` provides different information criteria to choose the tuning parameter \eqn{\lambda}. `cv.coxtp` uses  the cross validation to choose the tuning parameter.
#'
#' @seealso \code{\link{IC}}, \code{\link{plot}}, \code{\link{get.tvcoef}} and \code{\link{baseline}}.
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
#' Gray, R. J. (1992) Flexible methods for analyzing survival data using splines, with applications to breast cancer prognosis.
#' \emph{Journal of the American Statistical Association}, \strong{87}: 942-951. 
#' \cr
#' 
#' Gray, R. J. (1994) Spline-based tests in survival analysis.
#' \emph{Biometrics}, \strong{50}: 640-652.
#' \cr
#' 
#' Luo, L., He, K. Wu, W., and Taylor, J. M., (2023) Using information criteria to select smoothing parameters when analyzing survival data with time-varying coefficient hazard models.
#' \cr
#' 
#' Wu, W., Taylor, J. M., Brouwer, A. F., Luo, L., Kang, J., Jiang, H., and He, K. (2022) Scalable proximal methods for cause-specific hazard modeling with time-varying coefficients.
#' \emph{Lifetime Data Analysis}, \strong{28(2)}: 194-218.
#' \cr
#' 
#' Wood, S. N. (2017) P-splines with derivative based penalties and tensor product smoothing of unevenly distributed data.
#' \emph{Statistics and Computing}, \strong{27(4)}: 985-989.
#' \cr
#' 
#' Perperoglou, A., le Cessie, S., and van Houwelingen, H. C. (2006) A fast routine for fitting Cox models with time varying effects of the covariates.
#' \emph{Computer Methods and Programs in Biomedicine}, \strong{81(2)}: 154-161.
#' \cr
#' 
#' 
#' 
coxtp <- function(event , z , time ,strata=NULL ,penalty="Smooth-spline", nsplines=8, 
                  lambda = c(0.1,1,10), degree=3L,
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
  
  
  res <- vector(mode = "list", length = length(lambda))
  for(lambda_i in c(1:length(lambda))) {
    res[[lambda_i]] <- coxtp.base(fmla, data_NR, nsplines=nsplines, spline=spline, ties=ties, stop=stop,
                               method = method, btr = btr,
                               lambda_spline = lambda[lambda_i], TIC_prox = TIC_prox, ord = ord, degree = degree,
                               tol = tol, iter.max = iter.max, tau= tau, parallel = parallel, threads = threads,
                               fixedstep = fixedstep,
                               ICLastOnly = InfoCrit)
    
    names(res)[lambda_i] <- paste0("lambda",lambda_i)
    
    cat(paste("lambda", lambda[lambda_i], "is done.\n"))
    
  }
  
  res$lambda.list <- lambda
  attr(res, "class") <- c("list", "coxtp")
  attr(res, "fmla") <- fmla
  attr(res, "data_NR") <- data_NR
  attr(res, "nsplines") <- nsplines
  attr(res, "spline") <- spline
  attr(res, "ties") <- ties
  attr(res, "stop") <- stop
  attr(res, "method") <- method
  attr(res, "btr") <- btr
  attr(res, "TIC_prox") <- TIC_prox
  attr(res, "ord") <- ord
  attr(res, "degree") <- degree
  attr(res, "tol") <- tol
  attr(res, "iter.max") <- iter.max
  attr(res, "tau") <- tau
  attr(res, "parallel") <- parallel
  attr(res, "threads")  <- threads
  attr(res, "tau") <- tau
  attr(res, "fixedstep") <- fixedstep
  attr(res, "InfoCrit") <- InfoCrit
  
  
  return(res)
  
  
}




