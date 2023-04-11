#' calculating information criteria from a `coxtp` object
#' 
#' This function is to calculate information criteria from a `coxtp` object to select the penalization tuning parameter.
#'
#' @param fit model from `coxtp`.
#' @param IC.prox when calculating information criteria, there might be numerical issues (e.g. the Hessian matrix is close to be singular).
#' In such cases warnings will be given. 
#' If `IC.prox = true`, we modified the diagonal of the Hessian matrix a bit, which can lead to more stable estimates.
#' Default is `FALSE`.
#' 
#' 
#' @return 
#' \item{model.mAIC}{an object with S3 class \code{"coxtp"} using mAIC to select the tuning parameter.}
#' \item{model.TIC}{an object with S3 class \code{"coxtp"} using TIC to select the tuning parameter.}
#' \item{model.GIC}{an object with S3 class \code{"coxtp"} using GIC to select the tuning parameter.}
#' \item{mAIC}{a sequence of mAIC values for the different tuning parameters `lambda` from \code{"coxtp"}.}
#' \item{TIC}{a sequence of TIC values for the different tuning parameters `lambda` from \code{"coxtp"}.}
#' \item{GIC}{a sequence of GIC values for the different tuning parameters `lambda` from \code{"coxtp"}.}
#' 
#' @export
#' 
#' @examples
#' data(ExampleData)
#' z <- ExampleData$x
#' time <- ExampleData$time
#' event <- ExampleData$event
#' fit <- coxtp(event = event, z = z, time = time)
#' IC  <- IC(fit)
#' 
#' 
#' @details 
#' In order to select the proper smoothing parameter, we utilize the idea of information criteria. 
#' We provide four different information criteria to select the optimal smoothing parameter \eqn{\lambda}.
#' Generally, mAIC, TIC and GIC select similar parameters and the difference of resulting estimates are barely noticeable.
#' See details in Luo et al. (2022).
#' 
#' 
IC <- function(fit, IC.prox, ...){
  
  if (class(fit)[1]!= "list" & class(fit)[2]!= "coxtp" ) stop("fit is not an output from coxtp")
  
  lambda_list = fit$lambda.list
  
  AIC_all <- NULL
  TIC_all <- NULL
  GIC_all <- NULL
  
  for (index_lambda in c(1:length(lambda_list))) {
    lambda_i = lambda_list[index_lambda]
    
    fit.tmp = fit[[index_lambda]]
    theta.tmp = fit.tmp$ctrl.pts
    data = attr(fit.tmp, "data")
    term.event   = attr(fit.tmp, "term.event")
    term.tv      = attr(fit.tmp, "term.tv")
    term.time    = attr(fit.tmp, "term.time") 
    SmoothMatrix = attr(fit.tmp, "SmoothMatrix")
    SplineType   = fit.tmp$SplineType
    control      = attr(fit.tmp, "control")
    count.strata = attr(fit.tmp, "count.strata")
    bases        = fit.tmp$bases
    ties         = attr(fit.tmp, "ties")
    
    if(ties != "Breslow"){
      IC <-  ICcpp(event = data[,term.event], Z_tv = as.matrix(data[,term.tv]), B_spline = as.matrix(bases), 
                   count_strata = count.strata,
                   theta = theta.tmp, 
                   lambda_i = lambda_i,
                   SmoothMatrix  = SmoothMatrix,
                   SplineType    = SplineType,
                   method=control$method, 
                   lambda=control$lambda,
                   factor=control$factor,
                   parallel=control$parallel, threads=control$threads)
    } else{
      IC <-  ICcpp_bresties(event = data[,term.event], time = data[,term.time],
                   Z_tv = as.matrix(data[,term.tv]), B_spline = as.matrix(bases), 
                   count_strata = count.strata,
                   theta = theta.tmp, 
                   lambda_i = lambda_i,
                   SmoothMatrix  = SmoothMatrix,
                   SplineType    = SplineType,
                   method=control$method, 
                   lambda=control$lambda,
                   factor=control$factor,
                   parallel=control$parallel, threads=control$threads)
    }

    
    AIC_all <- c(AIC_all, IC$AIC)
    TIC_all <- c(TIC_all, IC$TIC)
    GIC_all <- c(GIC_all, IC$GIC)
  }
  

  AIC_lambda <- which.min(AIC_all)
  TIC_lambda <- which.min(TIC_all)
  GIC_lambda <- which.min(GIC_all)

  model_AIC <- fit[[AIC_lambda]]
  model_TIC <- fit[[TIC_lambda]]
  model_GIC <- fit[[GIC_lambda]]
  
  
  res <- NULL
  res$model.mAIC <- model_AIC
  res$model.TIC <- model_TIC
  res$model.GIC <- model_GIC
  res$mAIC <- AIC_all
  res$TIC <- TIC_all
  res$GIC <- GIC_all
  
  return(res)
  
  
  
}


