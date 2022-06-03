
#' Title
#'
#' Description
#'
#' @param event indicator for occurrence of event (1 if event, 0 else)
#' @param z define...
#' @param time = time to event
#' @param lambda_spline define...= lambda_spline,
#' @param spline define...= "Smooth-spline", 
#' @param nsplines define...= 8, 
#' @param ties define...= "none", 
#' @param tol define...= 1e-6, 
#' @param iter.max define...= 20L, 
#' @param method define...= "Newton",
#' @param btr define...= "dynamic", 
#' @param stop define...= "ratch", 
#' @param parallel define...= TRUE, 
#' @param threads define...= 3L, 
#' @param degree define...= 3L,
#' @param ICLastOnly define...= TRUE
#' 
#' 
#' @return 
#' A list containing ....
#' 
#' @examples
#' 
#' # load example data
#' data(simulN2kOP2)
#' 
#' z = simulN2kOP2[, c("z.1", "z.2")]
#' delta = simulN2kOP2$delta
#' time = simulN2kOP2$time
#' 
#' # specify smoothing parameter lambda
#' lambda_spline = c(1:5)
#' 
#' # model fitting
#' models = surtvep(event = delta, z = z, time = time, 
#'                   lambda_spline = lambda_spline,
#'                   spline="Smooth-spline", nsplines=8, ties="none", 
#'                   tol=1e-6, iter.max=20L, method="Newton",
#'                   btr="dynamic", stop="ratch", 
#'                   parallel=TRUE, threads=3L, degree=3L,
#'                   ICLastOnly = TRUE)
#' 
#' @export
surtvep <- function(event = delta, z = z, time = time, spline="Smooth-spline", nsplines=8, ties="Breslow", 
                    tol=1e-9, iter.max=20L, method="Newton", lambda=1e8,
                    btr="static", tau=0.5,
                    stop="ratch", parallel=FALSE, threads=1L, degree=3L, TIC = FALSE, TIC_prox = FALSE, 
                    lambda_spline = 0, ord = 4, fixedstep = FALSE, 
                    ICLastOnly = TRUE){
  
  
  stratum=rep(1, length(time))
  data_NR <- data.frame(event=delta, time=time, z, strata=stratum, stringsAsFactors=F)
  #Z.char <- paste0("X", 1:p)
  Z.char <- colnames(data_NR)[c(3:(ncol(data_NR)-1))]
  fmla <- formula(paste0("Surv(time, event)~",
                         paste(c(paste0("tv(", Z.char, ")"), "strata(strata)"), collapse="+")))
  
  # 
  
  lambda_all <- lambda_spline
  model <- list()
  
  for(lambda_index in 1:length(lambda_all)){
    
    model1 <- surtiver(fmla, data_NR, nsplines=10, spline=spline, ties=ties, stop=stop,
                       method = method, btr = btr,
                       lambda_spline = lambda_all[lambda_index],TIC_prox = TIC_prox, ord = ord, degree = degree, 
                       tol = tol, iter.max = iter.max, tau= tau, parallel = parallel, threads = threads,
                       fixedstep = FALSE,
                       penalizestop = FALSE,
                       ICLastOnly = ICLastOnly)
    
    model[[lambda_index]] <- model1
  }
  
  
  
  return(models = model)
}
