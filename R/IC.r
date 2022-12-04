#' calculating information criteria from a `coxtp` object
#' 
#' This function is to calculate information criteria from a `coxtp` object to select the penalization tuning paramter.
#'
#' @param fit model from `coxtp`.
#' @param IC.prox when calculating information criteria, there might be numerical issue(second order derivative, 
#' Hessian matrix approximate is singular). In such cases warnings will be given. 
#' if `true`, we modified the diagonal of hessian matrix a bit, which may lead to bias issues. 
#' Default is `FALSE`.
#' 
#' 
#' @return 
#' \item{model.AIC}{an object with S3 class \code{"coxtp"} using AIC to select the tunning parameter.}
#' \item{model.TIC}{an object with S3 class \code{"coxtp"} using TIC to select the tunning parameter.}
#' \item{model.GIC}{an object with S3 class \code{"coxtp"} using GIC to select the tunning parameter.}
#' \item{AIC}{a sequence of AIC values for the different tuning parameters `lambda` from \code{"coxtp"}.}
#' \item{TIC}{a sequence of AIC values for the different tuning parameters `lambda` from \code{"coxtp"}.}
#' \item{GIC}{a sequence of AIC values for the different tuning parameters `lambda` from \code{"coxtp"}.}
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
#' Generally, AIC, TIC and GIC selects similar parameters and the difference of resulting estimations are barely noticable.
#' See details in Lingfeng Luo et al. (2022)
#' 
#' 
IC <- function(fit, IC.prox, ...){
  
  if (class(fit.smoothspline)[1]!= "list" & class(fit.smoothspline)[2]!= "coxtp" ) stop("fit is not an output from coxtp")
  
  lambda_all = fit.smoothspline$lambda.list
  
  fmla <- attr(fit, "fmla") 
  data_NR<- attr(fit, "data_NR") 
  nsplines<- attr(fit, "nsplines") 
  spline<- attr(fit, "spline") 
  method<- attr(fit, "method") 
  btr<- attr(fit, "btr") 
  TIC_prox<- attr(fit, "TIC_prox") 
  ord<- attr(fit, "ord")
  degree<- attr(fit, "degree") 
  tol<- attr(fit, "tol") 
  iter.max<- attr(fit, "iter.max") 
  tau<- attr(fit, "tau")
  parallel<- attr(fit, "parallel") 
  fixedstep<- attr(fit, "fixedstep")
  InfoCrit<- attr(fit, "InfoCrit")
  ties<- attr(fit, "ties")
  stop<- attr(fit, "stop")
  threads <-  attr(fit, "threads")
  
  model <- vector(mode = "list", length = length(lambda_all))
  for(lambda_index in 1:length(lambda_all)){

    model1 <- coxtp.base(fmla, data_NR, nsplines=nsplines, spline=spline, ties=ties, stop=stop,
                         method = method, btr = btr,
                         lambda_spline = lambda_all[lambda_index],TIC_prox = TIC_prox, ord = ord, degree = degree,
                         tol = tol, iter.max = iter.max, tau= tau, parallel = parallel, threads = threads,
                         fixedstep = fixedstep,
                         ICLastOnly = TRUE)

    model[[lambda_index]] <- model1
  }


  AIC_all <- NULL
  TIC_all <- NULL
  GIC_all <- NULL

  for (i in 1:length(lambda_all)){
    AIC_all<-c(AIC_all, model[[i]]$AIC_all)
    TIC_all<-c(TIC_all, model[[i]]$TIC_all)
    GIC_all<-c(GIC_all, model[[i]]$GIC_all)
  }

  AIC_lambda <- which.min(AIC_all)
  TIC_lambda <- which.min(TIC_all)
  GIC_lambda <- which.min(GIC_all)

  model_AIC <- model[[AIC_lambda]]
  model_TIC <- model[[TIC_lambda]]
  model_GIC <- model[[GIC_lambda]]
  
  
  res <- NULL
  res$model.AIC <- model_AIC
  res$model.TIC <- model_TIC
  res$model.GIC <- model_GIC
  res$AIC <- AIC_all
  res$TIC <- TIC_all
  res$GIC <- GIC_all
  
  return(res)
  
  
  
}


