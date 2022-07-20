#' Example data for surtvep with 5000 observations of 4 variables:
#'
#' @format A data.frame 
#' \describe{
#'   \item{V1}{Simulated covariate V1, Binary variable with 0, 1. True time-dependent function is b(t)=1}
#'   \item{V2}{Simulated covariate V2, Binary variable with 0, 1. True time-dependent function is b(t)=exp(-1.5*t)}
#'   \item{event}{Simulated event variable, Binary varible with 0, 1 }
#'   \item{time}{Simulated time variable, Continous variable with non-negative value}
#' }
"sim_data"