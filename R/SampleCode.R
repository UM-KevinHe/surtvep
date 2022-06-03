
# rm(list=ls())  #remove list
# 
# library(RcppArmadillo)
# library(Rcpp)
# library(mvtnorm)
# library(splines)
# library(survival)
# library(ggplot2)
# library(ggrepel)

# load("simulN2kOP2.RData")
# 
# Rcpp::sourceCpp("PenalizeStopCpp.cpp")
# source("PenalizeStop.R")
# source("surtvp.R")
# 
# 
# #specify smoothing parameter lambda:
# lambda_spline = c(1:5)
# 
# ##model fitting:
# models <- surtvep(event = delta, z = z, time = time, 
#                   lambda_spline = lambda_spline,
#                   spline="Smooth-spline", nsplines=8, ties="none", 
#                   tol=1e-6, iter.max=20L, method="Newton",
#                   btr="dynamic", stop="ratch", 
#                   parallel=TRUE, threads=3L, degree=3L,
#                   ICLastOnly = TRUE)

