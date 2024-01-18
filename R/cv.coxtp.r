#' fit a cross-validated Cox non-proportional hazards model with P-spline or Smoothing-spline where penalization tuning parameter is provided by cross-validation
#' 
#' Fit a Cox non-proportional hazards model via penalized maximum likelihood. The penalization tuning parameter is provided by cross-validation.
#' 
#' @param event failure event response variable of length `nobs`, where `nobs` denotes the number of observations. It should be a vector containing 0 or 1.
#' @param z input covariate matrix, with `nobs` rows and `nvars` columns; each row is an observation. 
#' @param time observed event times, which should be a vector with non-negative values.
#' @param strata a vector of indicators for stratification. 
#' Default = `NULL` (i.e. no stratification group in the data), an unstratified model is implemented.
#' 
#' @param penalty a character string specifying the spline term for the penalized Newton method. 
#' This term is added to the log-partial likelihood, and the penalized log-partial likelihood serves as the new objective function to 
#' control the smoothness of the time-varying coefficients.
#' Default is `P-spline`. Three options are `P-spline`, `Smooth-spline` and `NULL`. 
#' If `NULL`, the method will be the same as `coxtv` (unpenalized time-varying effects models) and `lambda` (defined below)
#' will be set as 0. 
#' 
#' `P-spline` stands for Penalized B-spline. It combines the B-spline basis with a discrete quadratic penalty on the difference of basis coefficients between adjacent knots. 
#' When `lambda` goes to infinity, the time-varying effects are reduced to be constant. 
#' 
#' `Smooth-spline` refers to the Smoothing-spline, the derivative-based penalties combined with B-splines. See `degree` for different choices.
#' When `degree=3`, we use the cubic B-spline penalizing the second-order derivative, which reduces the time-varying effect to a linear term when `lambda` goes to infinity.
#' When `degree=2`, we use the quadratic B-spline penalizing first-order derivative, which reduces the time-varying effect to a constant when `lambda` goes to infinity. See Wood (2017) for details.
#' 
#' If `P-spline` or `Smooth-spline`, then `lambda` is initialized as a sequence (0.1, 1, 10). Users can modify `lambda`. See details in `lambda`.
#' 
#' @param lambda a user specified sequence as the penalization coefficients in front of the spline term specified by `penalty`. 
#' This is the tuning parameter for penalization.  The function `IC` can be used to select the best tuning parameter based on the information criteria. 
#' Users can specify larger values when the absolute values of the estimated time-varying effects are too large.
#' When `lambda` is `0`, Newton method without penalization is fitted. 
#' 
#' @param nfolds number of folds for cross-validation, the default value is 5. The smallest value allowable is `nfolds`=3.
#' @param foldid an optional vector of values between 1 and `nfolds` identifying what fold each observation is in. If supplied, `nfolds` can be missing.
#'
#' @param nsplines number of basis functions in the splines to span the time-varying effects. The default value is 8. 
#' We use the R function `splines::bs` to generate the B-splines. 
#' @param knots the internal knot locations (breakpoints) that define the B-splines.
#' The number of the internal knots should be `nsplines`-`degree`-1.
#' If `NULL`, the locations of knots are chosen as quantiles of distinct failure time points.
#' This choice leads to more stable results in most cases.
#' Users can specify the internal knot locations by themselves.
#' 
#' @param degree degree of the piecewise polynomial for generating the B-spline basis functions---default is 3 for cubic splines. 
#' `degree = 2` results in the quadratic B-spline basis functions. 
#' 
#' If `penalty` is `P-spline` or `NULL`, `degree`'s default value is 3. 
#' 
#' If `penalty` is `Smooth-spline`, `degree`'s default value is 2. 
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
#' @param iter.max maximum iteration number if the stopping criterion specified by `stop` is not satisfied. The default value is  `20`. 
#' @param method a character string specifying whether to use Newton method or proximal Newton method.  If `"Newton"` then Hessian is used, 
#' while the default method `"ProxN"` implements the proximal Newton which can be faster and more stable when there exists ill-conditioned second-order information of the log-partial likelihood.
#' See details in Wu et al. (2022).
#' 
#' @param gamma parameter for proximal Newton method `"ProxN"`. The default value is `1e8`.
#' @param btr a character string specifying the backtracking line-search approach. `"dynamic"` is a typical way to perform backtracking line-search. 
#' See details in Convex Optimization by Boyd and Vandenberghe (2004).
#' `"static"` limits Newton's increment and can achieve more stable results in some extreme cases, such as ill-conditioned second-order information of the log-partial likelihood, 
#' which usually occurs when some predictors are categorical with low frequency for some categories. 
#' Users should be careful with `static`, as this may lead to under-fitting.
#' @param tau a positive scalar used to control the step size inside the backtracking line-search. The default value is `0.5`.
#' @param parallel if `TRUE`, then the parallel computation is enabled. The number of threads in use is determined by `threads`.
#' @param threads an integer indicating the number of threads to be used for parallel computation. The default value is `2`. If `parallel` is false, then the value of `threads` has no effect.
#' @param fixedstep if `TRUE`, the algorithm will be forced to run `iter.max` steps regardless of the stopping criterion specified.
#' 
#' @return An object of class `"cv.coxtp"` is returned, which is a list with the ingredients of the cross-validation fit.
#' \item{model.cv}{a \code{"coxtp"} object with tuning parameter chosen based on cross-validation.} 
#' \item{lambda}{the values of `lambda` used in the fits.}
#' \item{cve}{the mean cross-validated error - a vector having the same length as lambda.
#' For the k-th testing fold (k = 1, ..., `nfolds`), we take the remaining folds as the training folds. 
#' Based on the model trained on the training folds, we calculate the log-partial likelihood on all the folds \eqn{loglik0} and training folds  \eqn{loglik1}. 
#' The `cve` is equal to \eqn{-2*(loglik0 - loglik1)}. See details in Verweij (1993). This approach avoids the construction of a partial likelihood on the test set so that the risk set is always sufficiently large.}
#' \item{lambda.min}{the value of `lambda` that gives minimum `cve`.}
#' 
#' @export
#' 
#'
#' @examples 
#' data(ExampleData)
#' z <- ExampleData$z
#' time  <- ExampleData$time
#' event <- ExampleData$event
#' lambda  = c(0.1, 1)
#' fit  <- cv.coxtp(event = event, z = z, time = time, lambda=lambda, nfolds = 5)
#' 
#' @details The function runs `coxtp` length of `lambda` by `nfolds` times; each is to compute the fit with each of the folds omitted.
#' 
#' @references 
#' Boyd, S., and Vandenberghe, L. (2004) Convex optimization. 
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
#' \emph{Statistical Methods in Medical Research}, \strong{in press}.
#' \cr
#' 
#' Perperoglou, A., le Cessie, S., and van Houwelingen, H. C. (2006) A fast routine for fitting Cox models with time varying effects of the covariates.
#' \emph{Computer Methods and Programs in Biomedicine}, \strong{81(2)}: 154-161.
#' \cr
#' 
#' Verweij, P. J., and Van Houwelingen, H. C. (1993) Cross‐validation in survival analysis.
#' \emph{Statistics in Medicine}, \strong{12(24)}: 2305-2314.
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
#' 
#'
cv.coxtp <- function(event , z, time, strata=NULL, 
                     lambda = c(0.1, 1, 10),
                     nfolds = 5, 
                     foldid = NULL, 
                     knots = NULL,
                     penalty="Smooth-spline", nsplines=8, ties="Breslow",
                     tol=1e-6, iter.max=20L, method="ProxN", gamma=1e8,
                     btr="dynamic", tau=0.5,
                     stop="ratch", parallel=FALSE, threads=1L, degree=3L, 
                     fixedstep = FALSE){
  
  if (nfolds <= 2) stop("nfolds must be greater or equal to 2!")
  
  spline = penalty
  ord = degree + 1
  TIC_prox = FALSE
  ICLastOnly = FALSE
  
  p = ncol(z)
  
  event  <- event[time_order <- order(time)]
  z      <- z[time_order,]
  strata <- strata[time_order]
  time   <- time[time_order]
  
  stratum <- if (length(strata) == 0) rep(1, length(time)) else strata
  
  #######################################################################
  ### 5-fold cross-validation
  #######################################################################
  N <- dim(z)[1]
  
  ind1 <- which(event==1)
  ind0 <- which(event==0)
  n1 <- length(ind1)
  n0 <- length(ind0)
  fold1 <- 1:n1 %% nfolds
  fold0 <- (n1 + 1:n0) %% nfolds
  fold1[fold1==0] <- nfolds
  fold0[fold0==0] <- nfolds
  fold <- integer(N)
  fold[event==1] <- sample(fold1)
  fold[event==0] <- sample(fold0)
  
  llk_difflambda_all <- NULL
  for (cv_fold in 1:nfolds) {
    cv_test_index = which(fold==cv_fold)
    
    delta_cv_train     = event[-cv_test_index]
    stratum_cv_train  = stratum[-cv_test_index]
    z_cv_train         = z[-cv_test_index,]
    time_cv_train      = time[-cv_test_index]
    knot_set = quantile(time_cv_train[delta_cv_train==1], prob=seq(1:(nsplines-4))/(nsplines-3))
    bs7      = splines::bs(time_cv_train,df=nsplines, knot=knot_set, intercept=TRUE, degree=degree)
    bs8      = matrix(bs7, nrow=dim(z_cv_train)[1]) # for build beta_t
    b_spline_cv_train = bs8
    
    delta_cv_test     = event
    stratum_cv_test  = stratum
    z_cv_test         = z
    time_cv_test      = time
    knot_set = quantile(time_cv_test[delta_cv_test==1], prob=seq(1:(nsplines-4))/(nsplines-3))
    bs7_test      = splines::bs(time_cv_test,df=nsplines, knot=knot_set, intercept=TRUE, degree=degree)
    bs8_test      = matrix(bs7_test, nrow=dim(z_cv_test)[1]) # for build beta_t
    b_spline_cv_test = bs8_test
    
    
    data_NR <- data.frame(event=delta_cv_train, time=time_cv_train, z_cv_train, strata=stratum_cv_train, stringsAsFactors=F)
    Z.char <- paste0("X", 1:p)
    fmla <- formula(paste0("Surv(time, event)~",
                           paste(c(paste0("tv(", Z.char, ")"), "strata(strata)"), collapse="+")))
    
    llk_difflambda <- NULL
    
    model <- list()
    for(lambda_index in 1:length(lambda)){
      
      model1 <- coxtp.base(fmla, data_NR, knots = knots, nsplines=nsplines, 
                           spline = spline, ties= ties, stop=stop,
                           TIC_prox = TIC_prox, ord = ord, degree = degree,
                           method = method, btr = btr, iter.max = 15, threads = threads, parallel = parallel,
                           lambda_spline = lambda[lambda_index], 
                           fixedstep = fixedstep)
      
      model[[lambda_index]] <- model1
      
      #model_all[[record_index]] <- model1
      llk_NR   <-  LogPartialTest(event = delta_cv_train, Z_tv =  z_cv_train, B_spline = as.matrix(b_spline_cv_train), 
                                  event_test = delta_cv_test, Z_tv_test = z_cv_test, B_spline_test = b_spline_cv_test,
                                  theta_list = model1$theta.list[length(model1$theta.list)],
                                  parallel = parallel, threads = threads, TestAll = TRUE)
      llk_difflambda <- c(llk_difflambda, -2*(length(delta_cv_test)*llk_NR$likelihood_all_test - length(delta_cv_train)*llk_NR$likelihood_all)) 
    }
    
    llk_difflambda_all <- rbind(llk_difflambda_all, llk_difflambda)
  }
  
  cve <- colMeans(llk_difflambda_all)
  
  lambda.min = lambda[which.min(cve)]
  
  data_NR <- data.frame(event=event, time=time, z, strata=stratum, stringsAsFactors=F)
  Z.char <- paste0("X", 1:p)
  fmla <- formula(paste0("Surv(time, event)~",
                         paste(c(paste0("tv(", Z.char, ")"), "strata(strata)"), collapse="+")))
  model.cv = coxtp.base(fmla, data_NR, knots = knots, nsplines=nsplines, spline = spline, ties= ties, stop=stop,
                        TIC_prox = TIC_prox, ord = ord, degree = degree,
                        method = method, btr = btr, iter.max = 15, threads = threads, parallel = parallel,
                        lambda_spline = lambda.min, 
                        fixedstep = fixedstep)
  
  res <- NULL
  res$cve = cve
  res$model.cv = model.cv
  res$lambda = lambda
  res$lambda.min= lambda.min
  
  class(res) <- "cv.coxtp"
  
  return(res)
  
  
}
