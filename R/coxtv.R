#' fit a Cox non-proportional hazards model
#' 
#' Fit a Cox non-proportional hazards model via maximum likelihood. 
#' 
#' @param event failure event response variable of length `nobs`, where `nobs` denotes the number of observations. It should be a vector containing 0 or 1.
#' @param z input covariate matrix, with `nobs` rows and `nvars` columns; each row is an observation. 
#' @param time observed event times, which should be a vector with non-negative values.
#' @param strata a vector of indicators for stratification. 
#' Default = `NULL` (i.e. no stratification group in the data), an unstratified model is implemented.
#' 
#' @param nsplines number of basis functions in the splines to span the time-varying effects. The default value is 8. 
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
#' @param ties a character string specifying the method for tie handling. If there are no tied events, 
#' the methods are equivalent.  
#' By default `"Breslow"` uses the Breslow approximation, which can be faster when many ties are present.
#' If `ties = "none"`, no approximation will be used to handle ties.
#' 
#' @param stop a character string specifying the stopping rule to determine convergence.  
#' `"incre"` means we stop the algorithm when Newton's increment is less than the `tol`. See details in Convex Optimization (Chapter 10) by Boyd and Vandenberghe (2004).
#' `"relch"` means we stop the algorithm when the \eqn{(loglik(m)-loglik(m-1))/(loglik(m))} is less than the `tol`,
#'  where \eqn{loglik(m)} denotes the log-partial likelihood at iteration step m.
#' `"ratch"` means we stop the algorithm when \eqn{(loglik(m)-loglik(m-1))/(loglik(m)-loglik(0))} is less than the `tol`.
#' `"all"` means we stop the algorithm when all the stopping rules (`"incre"`, `"relch"`, `"ratch"`) are met. 
#' The default value is `ratch`. 
#' If `iter.max` is achieved, it overrides any stop rule for algorithm termination.
#' 
#' @param tol tolerance used for stopping the algorithm. See details in `stop` below.
#'  The default value is  `1e-6`.
#' @param iter.max maximum iteration number if the stopping criterion specified by `stop` is not satisfied. The default value is  20.
#' @param method a character string specifying whether to use Newton method or proximal Newton method.  If `"Newton"` then Hessian is used, 
#' while the default method `"ProxN"` implements the proximal Newton which can be faster and more stable when there exists ill-conditioned second-order information of the log-partial likelihood.
#' See details in Wu et al. (2022).
#' @param gamma parameter for proximal Newton method `"ProxN"`. The default value is `1e8`.
#' @param btr a character string specifying the backtracking line-search approach. `"dynamic"` is a typical way to perform backtracking line-search. 
#' See details in Convex Optimization by Boyd and Vandenberghe (2004).
#' `"static"` limits Newton's increment and can achieve more stable results in some extreme cases, such as ill-conditioned second-order information of the log-partial likelihood, 
#' which usually occurs when some predictors are categorical with low frequency for some categories. 
#' Users should be careful with `static`, as this may lead to under-fitting.
#' @param tau a positive scalar used to control the step size inside the backtracking line-search. The default value is 0.5.
#' @param parallel if `TRUE`, then the parallel computation is enabled. The number of threads in use is determined by `threads`.
#' @param threads an integer indicating the number of threads to be used for parallel computation. The default value is `2`. If `parallel` is false, then the value of `threads` has no effect.
#' @param fixedstep if `TRUE`, the algorithm will be forced to run `iter.max` steps regardless of the stopping criterion specified.
#'
#'
#'
#' @return An object with S3 class `coxtv`.
#' \item{call}{the call that produced this object.}
#' \item{beta}{the estimated time-varying coefficient for each predictor at each unique time. 
#' It is a matrix of dimension `len_unique_t` by `nvars`, 
#' where `len_unique_t` is the length of unique observed event `time`s.}
#' 
#' \item{bases}{the basis matrix used in model fitting. If `ties="none"`, the dimension of the basis matrix is `nvars` by `nsplines`; 
#' if `ties="Breslow"`, the dimension is `len_unique_t` by `nsplines`. The matrix is constructed using the `bs::splines` function.}
#' \item{ctrl.pts}{estimated coefficient of the basis matrix of dimension `nvars` by `nsplines`. 
#' Each row represents a covariate's coefficient on the `nsplines`-dimensional basis functions.}
#' \item{Hessian}{the Hessian matrix of the log-partial likelihood, of which the dimension is `nsplines * nvars` by `nsplines * nvars`.}
#' \item{internal.knots}{the internal knot locations (breakpoints) that define the B-splines.}
#' \item{nobs}{number of observations.}
#' \item{theta.list}{the history of `ctrl.pts` of length `m` (the length of algorithm iterations), including `ctrl.pts` for each algorithm iteration.}
#' \item{VarianceMatrix}{the variance matrix of the estimated coefficients of the basis matrix, 
#' which is the inverse of the negative Hessian matrix.}
#' 
#' @export 
#' 
#' @seealso \code{\link{coef}}, \code{\link{plot}}, and the \code{\link{coxtp}} function.
#'
#' @examples 
#' data(ExampleData)
#' z <- ExampleData$z
#' time <- ExampleData$time
#' event <- ExampleData$event
#' fit <- coxtv(event = event, z = z, time = time)
#' 
#' @useDynLib surtvep, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats qnorm formula quantile
#' 
#' 
#' 
#' @details 
#' The model is fit by Newton method (proximal Newton method).
#' 
#' @references 
#' Boyd, S., Vandenberghe, L. (2004) Convex optimization. 
#' \emph{Cambridge University Press}.
#' \cr
#' 
#' Gray, R. J. (1992) Flexible methods for analyzing survival data using splines, with applications to breast cancer prognosis.
#' \emph{Journal of the American Statistical Association}, \strong{87(420)}: 942-951. 
#' \cr
#' 
#' Gray, R. J. (1994) Spline-based tests in survival analysis.
#' \emph{Biometrics}, \strong{50(3)}: 640-652.
#' \cr
#' 
#' Luo, L., He, K., Wu, W., and Taylor, J. M. (2023) Using information criteria to select smoothing parameters when analyzing survival data with time-varying coefficient hazard models.
#' \cr
#' 
#' Perperoglou, A., le Cessie, S., and van Houwelingen, H. C. (2006) A fast routine for fitting Cox models with time varying effects of the covariates.
#' \emph{Computer Methods and Programs in Biomedicine}, \strong{81(2)}: 154-161.
#' \cr
#' 
#' Wu, W., Taylor, J. M., Brouwer, A. F., Luo, L., Kang, J., Jiang, H., and He, K. (2022) Scalable proximal methods for cause-specific hazard modeling with time-varying coefficients.
#' \emph{Lifetime Data Analysis}, \strong{28(2)}: 194-218.
#' \cr
#' 
#'
#' 

coxtv <- function(event , z , time ,strata=NULL,  nsplines=8, 
                  knots = NULL, degree = 3,
                  ties="Breslow",
                  stop = 'ratch', tol = 1e-6, iter.max = 20, method = "ProxN", 
                  gamma = 1e8, btr = "dynamic",
                  tau = 0.5,
                  parallel=FALSE, threads=2L,fixedstep = FALSE){
  
  # order the data by time
  event  <- event[time_order <- order(time)]
  z      <- z[time_order,]
  strata <- strata[time_order]
  time   <- time[time_order]
  
  stratum <- if (length(strata) == 0) rep(1, length(time)) else strata
  
  data_NR <- data.frame(event=event, time=time, z, strata=stratum, stringsAsFactors=F)
  Z.char <- colnames(data_NR)[c(3:(ncol(data_NR)-1))]
  fmla <- formula(paste0("Surv(time, event)~",
                         paste(c(paste0("tv(", Z.char, ")"), "strata(strata)"), collapse="+")))

  spline="P-spline" # lambda = 0 so it doesn't matter
  
  res <- coxtv.base(fmla, data_NR, 
                    spline=spline,
                    nsplines = nsplines, ties=ties, stop=stop,
                    knots=knots,
                    method = method, btr = btr,
                    lambda_spline = 0, degree = degree,
                    tol = tol, iter.max = iter.max, tau= tau, parallel = parallel, threads = threads,
                    fixedstep = fixedstep)
    

  attr(res, "fmla") <- fmla
  attr(res, "data_NR") <- data_NR
  attr(res, "nsplines") <- nsplines
  attr(res, "ties") <- ties
  attr(res, "stop") <- stop
  attr(res, "method") <- method
  attr(res, "btr") <- btr
  attr(res, "degree") <- degree
  attr(res, "tol") <- tol
  attr(res, "iter.max") <- iter.max
  attr(res, "tau") <- tau
  attr(res, "parallel") <- parallel
  attr(res, "threads")  <- threads
  attr(res, "tau") <- tau
  attr(res, "fixedstep") <- fixedstep

  return(res)
}




