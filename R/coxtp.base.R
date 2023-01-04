coxtp.base <- function(formula, data, spline="Smooth-spline", nsplines=8, ties="Breslow", 
                     control, ...) {
  if (!ties%in%c("Breslow", "none")) stop("Invalid ties!")
  # pass ... args to coxtp.control
  extraArgs <- list(...)
  if (length(extraArgs)) {
    controlargs <- names(formals(coxtp.control))
    indx <- pmatch(names(extraArgs), controlargs, nomatch=0L)
    if (any(indx==0L))
      stop(gettextf("Argument(s) %s not matched!", 
                    names(extraArgs)[indx==0L]), domain=NA)
  }
  if (missing(control)) control <- coxtp.control(...)
  degree <- control$degree
  
  TIC           <- control$TIC
  TIC_prox      <- control$TIC_prox
  penalize      <- control$penalize
  lambda_spline <- control$lambda_spline
  ord           <- control$ord
  fixedstep     <- control$fixedstep
  addsecond     <- control$addsecond
  ICLastOnly    <- control$ICLastOnly
  
  
  
   # formula=fmla
  
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
  if (!spline%in%c("P-spline", "Smooth-spline") | 
      is.na(suppressWarnings(as.integer(nsplines[1]))) |
      as.integer(nsplines[1])<=degree+1) 
    stop(sprintf("Invalid spline or nsplines (should be at least %.0f)!", 
                 degree+2))
  nsplines <- nsplines[1]
  # model fitting
  if (spline=="P-spline") {
    knots <- 
      quantile(data[data[,term.event]==1,term.time], 
               (1:(nsplines-degree-1))/(nsplines-degree))
    
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
                                                  lambda_spline = lambda_spline,
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
                                                  ICLastOnly = control$ICLastOnly)
      # row.names(fit$ctrl.pts) <- term.tv
      # fit$internal.knots <- unname(knots)
      fit$uniqfailtimes <- uniqfailtimes.str
      fit$bases <- bases
      # return(fit)
    } else if (ties=="none") {
      bases <- 
        splines::bs(data[,term.time], degree=degree, intercept=T, 
                    knots=knots, Boundary.knots=range(times))
      fit <- 
        surtiver_fixtra_fit_penalizestop(data[,term.event], count.strata, 
                                         as.matrix(data[,term.tv]), as.matrix(bases), 
                                         matrix(0, length(term.tv), nsplines), 
                                         matrix(0, 1), matrix(0, 1),
                                         lambda_spline = lambda_spline,
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
                                         ICLastOnly = control$ICLastOnly)
      
      
      # row.names(fit$ctrl.pts) <- term.tv
      # fit$internal.knots <- unname(knots)
      fit$uniqfailtimes <- times
      fit$bases <- bases
      # return(fit)
    }
  } else if (spline=="Smooth-spline"){
    if (ties=="Breslow"){
      knots <- 
        quantile(data[data[,term.event]==1,term.time], 
                 (1:(nsplines-degree-1))/(nsplines-degree))
      uniqfailtimes.str <- unname(unlist(lapply(
        split(data[data[,term.event]==1,term.time],
              data[data[,term.event]==1,term.str]), unique)))
      bases <- splines::bs(uniqfailtimes.str, degree=degree, intercept=T, 
                           knots=knots, Boundary.knots=range(times))
      p_diffm   <- 1
      time      <-data[,term.time]
      x_seq     <- as.vector(c(min(time), knots, max(time)))
      h_j       <- diff(x_seq)
      #step 1:
      x_prime   <- x_seq
      if(ord == 4){
        knot_set2 <- c(min(time)-1,min(time)-1,min(time)-1, x_prime, max(time)+1, max(time)+1,max(time)+1)
        G_matrix  <- splines::splineDesign(knots = knot_set2 , x=x_prime ,ord = 4, derivs = 2)
      } else if(ord == 3){
        knot_set2 <- c(min(time)-1,min(time)-1, x_prime, max(time)+1, max(time)+1)
        G_matrix  <- splines::splineDesign(knots = knot_set2 , x=x_prime ,ord = 3, derivs = 1)
      } else {
        stop("ord must be 3 or 4")
      }
      P_matrix  <- matrix(0,p_diffm+1,p_diffm+1)
      H_matrix  <- matrix(0,p_diffm+1,p_diffm+1)
      for (i in 1:(p_diffm+1)) {
        for (j in 1:(p_diffm+1)) {
          P_matrix[i,j] <- (-1 + 2*(i-1)/p_diffm)^j
          H_matrix[i,j] <- (1 + (-1)^(i+j-2))/(i+j-1)
        }
      }
      W_tilde   <- t(solve(P_matrix))%*%H_matrix%*%solve(P_matrix)
      W_matrix  <- matrix(0,dim(G_matrix)[1],dim(G_matrix)[1])
      W_q       <- matrix(0,dim(G_matrix)[1],dim(G_matrix)[1])
      for(q in 1:length(h_j)){
        for (i in 1:(p_diffm+1)) {
          for (j in 1:(p_diffm+1)) {
            W_matrix[i+p_diffm*q-p_diffm,j+p_diffm*q-p_diffm] = 
              W_matrix[i+p_diffm*q-p_diffm,j+p_diffm*q-p_diffm] + h_j[q]*W_tilde[i,j]/2
            #W_matrix  <- W_matrix + W_q
          }
        }
      }
      SmoothMatrix  <- t(G_matrix)%*%W_matrix%*%G_matrix
      fit <- 
        surtiver_fixtra_fit_penalizestop_bresties(data[,term.event], data[,term.time], 
                                                  count.strata, 
                                                  as.matrix(data[,term.tv]), as.matrix(bases), 
                                                  matrix(0, length(term.tv), nsplines), 
                                                  matrix(0, 1), matrix(0, 1),
                                                  lambda_spline = lambda_spline,
                                                  SmoothMatrix  = SmoothMatrix,
                                                  effectsize    = control$effectsize,
                                                  SplineType    = "smooth-spline",
                                                  method=control$method, 
                                                  lambda=control$lambda, factor=control$factor,
                                                  parallel=control$parallel, threads=control$threads, 
                                                  tol=control$tol, iter_max=control$iter.max, 
                                                  s=control$s, t=control$t, 
                                                  btr=control$btr, stop=control$stop, TIC_prox = TIC_prox,
                                                  fixedstep=control$fixedstep,
                                                  difflambda = control$difflambda,
                                                  ICLastOnly = control$ICLastOnly)
      # row.names(fit$ctrl.pts) <- term.tv
      # fit$internal.knots <- unname(knots)
      fit$uniqfailtimes <- uniqfailtimes.str
      fit$bases <- bases
      fit$knots <- knots
      # return(fit)
      
    } else if (ties == "none"){
      knots <- 
        quantile(data[data[,term.event]==1,term.time], 
                 (1:(nsplines-degree-1))/(nsplines-degree))
      bases <- 
        splines::bs(data[,term.time], degree=degree, intercept=T, 
                    knots=knots, Boundary.knots=range(times))
      p_diffm   <- 1
      time      <-data[,term.time]
      x_seq     <- as.vector(c(min(time), knots, max(time)))
      h_j       <- diff(x_seq)
      #step 1:
      x_prime   <- x_seq
      if(ord == 4){
        knot_set2 <- c(min(time)-1,min(time)-1,min(time)-1, x_prime, max(time)+1, max(time)+1,max(time)+1)
        G_matrix  <- splines::splineDesign(knots = knot_set2 , x=x_prime ,ord = 4, derivs = 2)
      } else if(ord == 3){
        knot_set2 <- c(min(time)-1,min(time)-1, x_prime, max(time)+1, max(time)+1)
        G_matrix  <- splines::splineDesign(knots = knot_set2 , x=x_prime ,ord = 3, derivs = 1)
      } else {
        stop("ord must be 3 or 4")
      }
      P_matrix  <- matrix(0,p_diffm+1,p_diffm+1)
      H_matrix  <- matrix(0,p_diffm+1,p_diffm+1)
      for (i in 1:(p_diffm+1)) {
        for (j in 1:(p_diffm+1)) {
          P_matrix[i,j] <- (-1 + 2*(i-1)/p_diffm)^j
          H_matrix[i,j] <- (1 + (-1)^(i+j-2))/(i+j-1)
        }
      }
      W_tilde   <- t(solve(P_matrix))%*%H_matrix%*%solve(P_matrix)
      W_matrix  <- matrix(0,dim(G_matrix)[1],dim(G_matrix)[1])
      W_q       <- matrix(0,dim(G_matrix)[1],dim(G_matrix)[1])
      for(q in 1:length(h_j)){
        for (i in 1:(p_diffm+1)) {
          for (j in 1:(p_diffm+1)) {
            W_matrix[i+p_diffm*q-p_diffm,j+p_diffm*q-p_diffm] = 
              W_matrix[i+p_diffm*q-p_diffm,j+p_diffm*q-p_diffm] + h_j[q]*W_tilde[i,j]/2
            #W_matrix  <- W_matrix + W_q
          }
        }
      }
      SmoothMatrix  <- t(G_matrix)%*%W_matrix%*%G_matrix
      fit <- 
        surtiver_fixtra_fit_penalizestop(data[,term.event], count.strata, 
                                         as.matrix(data[,term.tv]), as.matrix(bases), 
                                         matrix(0, length(term.tv), nsplines), 
                                         matrix(0, 1), matrix(0, 1),
                                         lambda_spline = lambda_spline,
                                         SmoothMatrix  = SmoothMatrix,
                                         effectsize    = control$effectsize,
                                         SplineType    = "smooth-spline",
                                         method=control$method, 
                                         lambda=control$lambda, factor=control$factor,
                                         parallel=control$parallel, threads=control$threads, 
                                         tol=control$tol, iter_max=control$iter.max, 
                                         s=control$s, t=control$t, 
                                         btr=control$btr, stop=control$stop, TIC_prox = TIC_prox,
                                         fixedstep=control$fixedstep,
                                         difflambda = control$difflambda,
                                         ICLastOnly = control$ICLastOnly)
      # row.names(fit$ctrl.pts) <- term.tv
      # fit$internal.knots <- unname(knots)
      fit$uniqfailtimes <- times
      fit$bases <- bases
      # return(fit)
    }
    
  }
  
  res <- NULL
  res$theta.list  <- fit$theta_list
  res$VarianceMatrix <- fit$VarianceMatrix
  res$times <- times
  res$ctrl.pts <- fit$theta_list[[length(fit$theta_list)]]
  res$internal.knots    <- unname(knots)
  res$uniqfailtimes     <- fit$uniqfailtimes.str
  res$bases <- fit$bases
  res$AIC_all   <- fit$AIC_all
  res$TIC_all   <- fit$TIC_all
  res$GIC_all   <- fit$GIC_all
  res$SplineType <- fit$SplineType
  res$info <- fit$info
  
  class(res) <- "coxtp"

  row.names(res$ctrl.pts)   <- term.tv
  attr(res, "data")         <- data
  attr(res, "term.event")   <- term.event
  attr(res, "term.tv")      <- term.tv
  attr(res, "term.time")    <- term.time
  attr(res, "SmoothMatrix") <- SmoothMatrix
  attr(res, "spline") <- spline
  attr(res, "count.strata") <- count.strata
  # if (length(term.ti)>0) {
  #   fit$tief <- c(fit$tief)
  #   names(fit$tief) <- term.ti
  # }
  # colnames(fit$info) <- rownames(fit$info) <-
  #   c(rep(term.tv, each=nsplines), term.ti)
  attr(res, "time")     <- time
  attr(res, "ties")     <- ties
  attr(res, "nsplines") <- nsplines
  attr(res, "degree.spline") <- degree
  attr(res, "control") <- control
  attr(res, "response") <- term.event
  attr(res, "internal.knots")    <- unname(knots)
  class(res) <- "coxtp"
  return(res)
}

coxtp.control <- function(tol=1e-9, iter.max=20L, method="Newton", lambda=1e8,
                             factor=10, btr="dynamic", sigma=1e-2, tau=0.6,
                             stop="incre", parallel=FALSE, threads=1L, degree=3L, TIC = FALSE, TIC_prox = FALSE, 
                             lambda_spline = 0, ord = 4, fixedstep = FALSE, effectsize = 0, difflambda = FALSE,
                             addsecond = FALSE,  ICLastOnly = FALSE) {
  if (tol <= 0) stop("Invalid convergence tolerance!")
  if (iter.max <= 0 | !is.numeric(iter.max))
    stop("Invalid maximum number of iterations!")
  if (method %in% c("Newton", "ProxN")) {
    if (method=="ProxN" & (lambda <= 0 | factor < 1))
      stop("Argument lambda <=0 or factor < 1 when method = 'ProxN'!")
  } else stop("Invalid estimation method!")
  
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
       parallel=parallel, threads=as.integer(threads), degree=as.integer(degree),
       TIC=TIC, TIC_prox = TIC_prox, lambda_spline= lambda_spline, ord = ord, fixedstep = fixedstep, 
       effectsize = effectsize, difflambda = difflambda, addsecond = addsecond, 
       ICLastOnly=ICLastOnly)
}


VarianceMatrix <- function(formula, data, spline="P-spline", nsplines=8, ties="Breslow", 
                           theta_opt_lambda, opt_lambda,
                           control,...){
  
  if (!ties%in%c("Breslow", "none")) stop("Invalid ties!")
  # pass ... args to coxtp.control
  extraArgs <- list(...)
  if (length(extraArgs)) {
    controlargs <- names(formals(coxtp.control))
    indx <- pmatch(names(extraArgs), controlargs, nomatch=0L)
    if (any(indx==0L))
      stop(gettextf("Argument(s) %s not matched!", 
                    names(extraArgs)[indx==0L]), domain=NA)
  }
  if (missing(control)) control <- coxtp.control(...)
  degree <- control$degree
  
  TIC           <- control$TIC
  TIC_prox      <- control$TIC_prox
  penalize      <- control$penalize
  lambda_spline <- control$lambda_spline
  ord           <- control$ord
  fixedstep     <- control$fixedstep
  addsecond     <- control$addsecond

  Terms <- terms(formula, specials=c("tv", "strata", "offset"))
  if(attr(Terms, 'response')==0) stop("Formula must have a Surv response!")
  factors <- attr(Terms, 'factors')
  terms <- row.names(factors)
  idx.r <- attr(Terms, 'response')
  idx.o <- attr(Terms, 'specials')$offset
  idx.str <- attr(Terms, 'specials')$strata
  idx.tv <- attr(Terms, 'specials')$tv
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
  if (!spline%in%c("P-spline", "Smooth-spline") | 
      is.na(suppressWarnings(as.integer(nsplines[1]))) |
      as.integer(nsplines[1])<=degree+1) 
    stop(sprintf("Invalid spline or nsplines (should be at least %.0f)!", 
                 degree+2))
  nsplines <- nsplines[1]
  # model fitting
  if (spline=="P-spline") {
    knots <- 
      quantile(data[data[,term.event]==1,term.time], 
               (1:(nsplines-degree-1))/(nsplines-degree))
    
    if (ties=="Breslow"){
      uniqfailtimes.str <- unname(unlist(lapply(
        split(data[data[,term.event]==1,term.time],
              data[data[,term.event]==1,term.str]), unique)))
      bases <- splines::bs(uniqfailtimes.str, degree=degree, intercept=T, 
                           knots=knots, Boundary.knots=range(times))
    } else if (ties=="none") {
      bases <- 
        splines::bs(data[,term.time], degree=degree, intercept=T, 
                    knots=knots, Boundary.knots=range(times))
    }
  } else if (spline=="Smooth-spline"){
    if (ties=="Breslow"){
      knots <- 
        quantile(data[data[,term.event]==1,term.time], 
                 (1:(nsplines-degree-1))/(nsplines-degree))
      uniqfailtimes.str <- unname(unlist(lapply(
        split(data[data[,term.event]==1,term.time],
              data[data[,term.event]==1,term.str]), unique)))
      bases <- splines::bs(uniqfailtimes.str, degree=degree, intercept=T, 
                           knots=knots, Boundary.knots=range(times))
      p_diffm   <- 1
      time      <-data[,term.time]
      x_seq     <- as.vector(c(min(time), knots, max(time)))
      h_j       <- diff(x_seq)
      #step 1:
      x_prime   <- x_seq
      if(ord == 4){
        knot_set2 <- c(min(time)-1,min(time)-1,min(time)-1, x_prime, max(time)+1, max(time)+1,max(time)+1)
        G_matrix  <- splines::splineDesign(knots = knot_set2 , x=x_prime ,ord = 4, derivs = 2)
      } else if(ord == 3){
        knot_set2 <- c(min(time)-1,min(time)-1, x_prime, max(time)+1, max(time)+1)
        G_matrix  <- splines::splineDesign(knots = knot_set2 , x=x_prime ,ord = 3, derivs = 1)
      } else {
        stop("ord must be 3 or 4")
      }
      P_matrix  <- matrix(0,p_diffm+1,p_diffm+1)
      H_matrix  <- matrix(0,p_diffm+1,p_diffm+1)
      for (i in 1:(p_diffm+1)) {
        for (j in 1:(p_diffm+1)) {
          P_matrix[i,j] <- (-1 + 2*(i-1)/p_diffm)^j
          H_matrix[i,j] <- (1 + (-1)^(i+j-2))/(i+j-1)
        }
      }
      W_tilde   <- t(solve(P_matrix))%*%H_matrix%*%solve(P_matrix)
      W_matrix  <- matrix(0,dim(G_matrix)[1],dim(G_matrix)[1])
      W_q       <- matrix(0,dim(G_matrix)[1],dim(G_matrix)[1])
      for(q in 1:length(h_j)){
        for (i in 1:(p_diffm+1)) {
          for (j in 1:(p_diffm+1)) {
            W_matrix[i+p_diffm*q-p_diffm,j+p_diffm*q-p_diffm] = 
              W_matrix[i+p_diffm*q-p_diffm,j+p_diffm*q-p_diffm] + h_j[q]*W_tilde[i,j]/2
            #W_matrix  <- W_matrix + W_q
          }
        }
      }
      SmoothMatrix  <- t(G_matrix)%*%W_matrix%*%G_matrix

    } else if (ties == "none"){
      knots <- 
        quantile(data[data[,term.event]==1,term.time], 
                 (1:(nsplines-degree-1))/(nsplines-degree))
      bases <- 
        splines::bs(data[,term.time], degree=degree, intercept=T, 
                    knots=knots, Boundary.knots=range(times))
      p_diffm   <- 1
      
      time      <-data[,term.time]
      x_seq     <- as.vector(c(min(time), knots, max(time)))
      h_j       <- diff(x_seq)
      #step 1:
      x_prime   <- x_seq
      if(ord == 4){
        knot_set2 <- c(min(time)-1,min(time)-1,min(time)-1, x_prime, max(time)+1, max(time)+1,max(time)+1)
        G_matrix  <- splines::splineDesign(knots = knot_set2 , x=x_prime ,ord = 4, derivs = 2)
      } else if(ord == 3){
        knot_set2 <- c(min(time)-1,min(time)-1, x_prime, max(time)+1, max(time)+1)
        G_matrix  <- splines::splineDesign(knots = knot_set2 , x=x_prime ,ord = 3, derivs = 1)
      } else {
        stop("ord must be 3 or 4")
      }
      P_matrix  <- matrix(0,p_diffm+1,p_diffm+1)
      H_matrix  <- matrix(0,p_diffm+1,p_diffm+1)
      for (i in 1:(p_diffm+1)) {
        for (j in 1:(p_diffm+1)) {
          P_matrix[i,j] <- (-1 + 2*(i-1)/p_diffm)^j
          H_matrix[i,j] <- (1 + (-1)^(i+j-2))/(i+j-1)
        }
      }
      W_tilde   <- t(solve(P_matrix))%*%H_matrix%*%solve(P_matrix)
      W_matrix  <- matrix(0,dim(G_matrix)[1],dim(G_matrix)[1])
      W_q       <- matrix(0,dim(G_matrix)[1],dim(G_matrix)[1])
      for(q in 1:length(h_j)){
        for (i in 1:(p_diffm+1)) {
          for (j in 1:(p_diffm+1)) {
            W_matrix[i+p_diffm*q-p_diffm,j+p_diffm*q-p_diffm] = 
              W_matrix[i+p_diffm*q-p_diffm,j+p_diffm*q-p_diffm] + h_j[q]*W_tilde[i,j]/2
            #W_matrix  <- W_matrix + W_q
          }
        }
      }
      SmoothMatrix  <- t(G_matrix)%*%W_matrix%*%G_matrix
      
      fit <-  VarianceMatrixCalculate(event = data[,term.event], Z_tv = as.matrix(data[,term.tv]), B_spline = as.matrix(bases), 
                                      theta = theta_opt_lambda, 
                                      lambda_i = opt_lambda,
                                      SmoothMatrix  = SmoothMatrix,
                                      SplineType    = "smooth-spline",
                                      method=control$method, 
                                      lambda=control$lambda,
                                      factor=control$factor,
                                      parallel=control$parallel, threads=control$threads)
      
      fit$bases <- bases
      return(fit)
    }
  
  
  
  }
}


#' @rdname confint.coxtv
#' 
#' @param fit fitted \code{"coxtp"} model
#' @param times the time interval to be estamtied. The default value is the time of the fitted model
#' @param parm the names of parameter
#' @param level the confidence level. Default is 0.95.
#' 
#' 
#' @examples 
#' \dontrun{
#' data(ExampleData)
#' z <- ExampleData$x
#' time <- ExampleData$time
#' event <- ExampleData$event
#' fit <- coxtp(event = event, z = z, time = time)
#' confint(fit$lambda1)
#' }
#' 
#' @exportS3Method confint coxtp
confint.coxtp <- function(fit, times, parm, level=0.95) {
  
  if (missing(fit)) stop ("Argument fit is required!")
  if (class(fit)!="coxtp") stop("Object fit is not of class 'coxtp'!")
  # if (missing(times)) {
  times <- fit$times
  # } else {
  #   if (!is.numeric(times) | min(times)<0) stop("Invalid times!")
  # }
  if (!is.numeric(level) | level[1]>1 | level[1]<0) stop("Invalid level!")
  level <- level[1]
  times <- times[order(times)]
  times <- unique(times)
  spline <- attr(fit, "spline"); degree <- attr(fit, "degree.spline")
  knots <- attr(fit, "internal.knots"); nsplines <- attr(fit, "nsplines")
  method <- attr(fit, "control")$method
  term.tv <- rownames(fit$ctrl.pts)
  # if (missing(parm)) {
    parm <- term.tv
  # } else if (length(parm)>0) {
  #   indx <- pmatch(parm, term.tv, nomatch=0L)
  #   if (any(indx==0L))
  #     stop(gettextf("%s not matched!", parm[indx==0L]), domain=NA)
  # } else stop("Invalid parm!")
  rownames.info <- rep(term.tv, each=nsplines)
  if (method=="Newton") {
    invinfo <- solve(fit$info)
    # invinfo <- fit$VarianceMatrix
  } else if (method=="ProxN") {
    invinfo <- solve(fit$info+diag(sqrt(.Machine$double.eps),dim(fit$VarianceMatrix)[1]))
    # invinfo <- fit$VarianceMatrix
  }
  # parm.ti <- intersect(parm, c(term.ti))
  parm.tv <- intersect(parm, c(term.tv))
  quant.upper <- qnorm((1+level)/2)
  ls <- list()
  # if (length(parm.ti)!=0) {
  #   est.ti <- fit$tief[term.ti%in%parm.ti]
  #   se.ti <- c(sqrt(diag(as.matrix(invinfo[rownames.info%in%parm.ti,
  #                                          rownames.info%in%parm.ti]))))
  #   mat.ti <- cbind(est.ti, est.ti-quant.upper*se.ti, est.ti+quant.upper*se.ti)
  #   colnames(mat.ti) <- 
  #     c("est", paste0(round(100*c(1-(1+level)/2,(1+level)/2),1),"%"))
  #   rownames(mat.ti) <- parm.ti
  #   ls$tief <- mat.ti
  # }
  # if (length(parm.tv)!=0) {
  # if (spline=="B-spline") {
  bases <- splines::bs(times, degree=degree, intercept=T, knots=knots, 
                       Boundary.knots=range(fit$times))
  ctrl.pts <- matrix(fit$ctrl.pts[term.tv%in%parm.tv,], ncol=nsplines)
  ls$tvef <- lapply(parm.tv, function(tv) {
    est.tv <- bases%*%ctrl.pts[parm.tv%in%tv,]
    se.tv <- sqrt(apply(bases, 1, function(r) {
      idx <- rownames.info%in%tv
      return(t(r)%*%invinfo[idx, idx]%*%r)}))
    mat.tv <- cbind(est.tv, est.tv-quant.upper*se.tv, 
                    est.tv+quant.upper*se.tv)
    colnames(mat.tv) <- 
      c("est", paste0(round(100*c(1-(1+level)/2,(1+level)/2),1),"%"))
    rownames(mat.tv) <- times
    return(mat.tv)
  })
  names(ls$tvef) <- parm.tv
  # } 
  # else if (spline=="P-spline") {
  
  # }
  # }
  return(ls)
}
