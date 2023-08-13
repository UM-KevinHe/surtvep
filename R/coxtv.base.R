#' @importFrom stats qnorm formula quantile

coxtv.base <- function(formula, data, spline, ties="Breslow", 
                  control, ...) {
  if (!ties%in%c("Breslow", "none")) stop("Invalid ties!")
  # pass ... args to coxtp.control
  extraArgs <- list(...)
  if (length(extraArgs)) {
    controlargs <- names(formals(coxtv.control))
    indx <- pmatch(names(extraArgs), controlargs, nomatch=0L)
    if (any(indx==0L))
      stop(gettextf("Argument(s) %s not matched!", 
                    names(extraArgs)[indx==0L]), domain=NA)
  }
  if (missing(control)) control <- coxtv.control(...)

  knots         <- control$knots
  degree        <- control$degree
  nsplines      <- control$nsplines
  TIC_prox      <- control$TIC_prox
  penalize      <- control$penalize
  lambda_spline <- control$lambda_spline
  ord           <- control$ord
  fixedstep     <- control$fixedstep

  Terms <- terms(formula, specials=c("tv", "strata", "offset"))
  if(attr(Terms, 'response')==0) stop("Formula must have a Surv response!")
  factors <- attr(Terms, 'factors')
  terms <- row.names(factors)
  idx.r <- attr(Terms, 'response')
  idx.o <- attr(Terms, 'specials')$offset
  idx.str <- attr(Terms, 'specials')$strata
  idx.tv <- attr(Terms, 'specials')$tv
  
  # #Add time:
  # time  <-data[,idx.str]
  
  if (is.null(idx.tv)) stop("No variable specified with time-variant effect!")
  term.ti <- terms[setdiff(1:length(terms), c(idx.r, idx.o, idx.str, idx.tv))]
  term.time <- gsub(".*\\(([^,]*),\\s+([^,]*)\\)", "\\1", terms[idx.r])
  term.event <- gsub(".*\\(([^,]*),\\s+([^,]*)\\)", "\\2", terms[idx.r])
  if (!is.null(idx.o)) term.o <- gsub(".*\\((.*)\\)", "\\1", terms[idx.o])
  term.str <- ifelse(!is.null(idx.str),
                     gsub(".*\\((.*)\\)", "\\1", terms[idx.str][1]), "strata")
  
  term.tv <- gsub(".*\\((.*)\\)", "\\1", terms[idx.tv])
  if (is.null(idx.str)) data[,term.str] <- "unstratified"
  data <- data[order(data[,term.str], data[,term.time], -data[,term.event]),]
  row.names(data) <- NULL
  times <- data[,term.time]; times <- unique(times); times <- times[order(times)]
  strata.noevent <- 
    unique(data[,term.str])[sapply(split(data[,term.event],data[,term.str]), 
                                   sum)==0]
  data <- data[!data[,term.str]%in%strata.noevent,] # drop strata with no event
  count.strata <- sapply(split(data[,term.str], data[,term.str]), length)
  
  
  if (any(!data[,term.event]%in%c(0,1)) |
      !is.numeric(data[,term.time]) |
      min(data[,term.time])<0) stop("Invalid Surv object!")
  # check spline-related arguments
  
  if (is.na(suppressWarnings(as.integer(nsplines[1]))) |
      as.integer(nsplines[1])<=degree+1) 
    stop(sprintf("Invalid nsplines (should be at least %.0f)!", 
                 degree+2))
  
  nsplines <- nsplines[1]
  
  
  if (spline=="P-spline") {
    if(is.null(knots)){
      knots <- 
        quantile(data[data[,term.event]==1,term.time], 
                 (1:(nsplines-degree-1))/(nsplines-degree))
    }
    if (ties=="Breslow"){
      uniqfailtimes.str <- unname(unlist(lapply(
        split(data[data[,term.event]==1,term.time],
              data[data[,term.event]==1,term.str]), unique)))
      bases <- splines::bs(uniqfailtimes.str, degree=degree, intercept=T, 
                           knots=knots, Boundary.knots=range(times))
      
      fit <- 
        surtiver_fixtra_fit_penalizestop_bresties(data[,term.event], data[,term.time], 
                                                  count.strata, 
                                                  as.matrix(data[,term.tv]), as.matrix(bases), 
                                                  matrix(0, length(term.tv), nsplines), 
                                                  matrix(0, 1), matrix(0, 1),
                                                  lambda_spline = 0,
                                                  SmoothMatrix  = matrix(0,1),
                                                  effectsize    = control$effectsize,
                                                  SplineType    = "pspline",
                                                  method=control$method, 
                                                  lambda=control$lambda, factor=control$factor,
                                                  parallel=control$parallel, threads=control$threads, 
                                                  tol=control$tol, iter_max=control$iter.max, 
                                                  s=control$s, t=control$t, 
                                                  btr=control$btr, stop=control$stop, TIC_prox = FALSE,
                                                  fixedstep = control$fixedstep,
                                                  difflambda = control$difflambda,
                                                  ICLastOnly = FALSE)
      fit$uniqfailtimes <- uniqfailtimes.str
      fit$bases <- bases
    } else if (ties=="none") {
      bases <- 
        splines::bs(data[,term.time], degree=degree, intercept=T, 
                    knots=knots, Boundary.knots=range(times))
      fit <- 
        surtiver_fixtra_fit_penalizestop(data[,term.event], count.strata, 
                                         as.matrix(data[,term.tv]), as.matrix(bases), 
                                         matrix(0, length(term.tv), nsplines), 
                                         matrix(0, 1), matrix(0, 1),
                                         lambda_spline = 0,
                                         SmoothMatrix  = matrix(0,1),
                                         effectsize    = control$effectsize,
                                         SplineType    = "pspline",
                                         method=control$method, 
                                         lambda=control$lambda, factor=control$factor,
                                         parallel=control$parallel, threads=control$threads, 
                                         tol=control$tol, iter_max=control$iter.max, 
                                         s=control$s, t=control$t, 
                                         btr=control$btr, stop=control$stop, TIC_prox = control$TIC_prox,
                                         fixedstep = control$fixedstep,
                                         difflambda = control$difflambda,
                                         ICLastOnly = FALSE)
      
      fit$uniqfailtimes <- times
      fit$bases <- bases
    }
  } else if (spline=="Smooth-spline"){
    
  }
  
  
  res <- NULL
  res$theta.list  <- fit$theta_list
  res$VarianceMatrix <- fit$VarianceMatrix
  res$times <- times
  res$ctrl.pts <- fit$theta_list[[length(fit$theta_list)]]
  res$internal.knots    <- unname(knots)
  res$uniqfailtimes     <- fit$uniqfailtimes.str
  res$bases <- fit$bases
  res$info <- fit$info
  
  row.names(res$ctrl.pts) <- term.tv
  
  
  class(res) <- "coxtv"
  attr(res, "spline")     <- spline
  attr(res, "count.strata")     <- count.strata
  attr(res, "basehazard") <- fit$hazard
  
  if (length(term.ti)>0) {
    res$tief <- c(res$tief)
    names(res$tief) <- term.ti
  }
  attr(res, "data") <- data

  attr(res, "nsplines") <- nsplines
  attr(res, "degree.spline") <- degree
  attr(res, "control") <- control
  attr(res, "response") <- term.event
  return(res)
}


coxtv.control <- function(tol=1e-9, iter.max=20L, method="ProxN", gamma=1e8,
                          knots=NULL, nsplines=8L,
                          factor=10, btr="dynamic", sigma=1e-2, tau=0.6,
                          stop="incre", parallel=FALSE, threads=1L, degree=3L, TIC_prox = FALSE, 
                          lambda_spline = 0, ord = 4, fixedstep = FALSE, difflambda = FALSE, effectsize = 0) {
  lambda = gamma
  
  if (tol <= 0) stop("Invalid convergence tolerance!")
  if (iter.max <= 0 | !is.numeric(iter.max))
    stop("Invalid maximum number of iterations!")
  if (method %in% c("Newton", "ProxN")) {
    if (method=="ProxN" & (lambda <= 0 | factor < 1))
      stop("Argument lambda <=0 or factor < 1 when method = 'ProxN'!")
  } else stop("Invalid estimation method!")
  if (!is.null(knots)){
    if(length(knots)!=nsplines-degree-1)
      stop("The number of the internal knots should be `nsplines`-`degree`-1!")
  }
  if (btr %in% c("none", "static", "dynamic")) {
    if (sigma <= 0 | tau <= 0) 
      stop("Search control parameter sigma <= 0 or sigma <= 0!")
    if (btr=="dynamic") {
      if (sigma > 0.5)
        stop("Search control parameter sigma > 0.5!")
      if (tau < 0.5) stop("Search control parameter tau < 0.5!")
    }
  } else stop("Invalid backtracking line search approach!")
  if (!is.logical(parallel))
    stop("Invalid argument parallel! It should be a logical constant.")
  if (!is.numeric(threads) | threads <= 0)
    stop("Invalid number of threads! It should be a positive integer.")
  if (parallel & threads==1L)
    stop(paste("Invalid number of threads for parallel computing!",
               "It should be greater than one."))
  if (!stop %in% c("incre", "ratch", "relch")) stop("Invalid stopping rule!")
  if (!is.numeric(degree) | degree < 2L) stop("Invalid degree for spline!")
  list(tol=tol, iter.max=as.integer(iter.max), method=method, lambda=lambda,
       factor=factor, btr=btr, s=sigma, t=tau, stop=stop,
       knots=knots,
       parallel=parallel, threads=as.integer(threads), degree=as.integer(degree),
       TIC_prox = TIC_prox, lambda_spline= lambda_spline, ord = ord, fixedstep = fixedstep, difflambda = difflambda,
       effectsize = effectsize, nsplines = nsplines)
}



