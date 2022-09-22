#' Example data for surtvep with 5000 observations of 8 variables(With strata):
#'
#' @format A data.frame 
#' \describe{
#' \itemize{
#'   \item{V1}{Simulated covariate V1, Binary variable with 0, 1. True time-dependent function is b(t)=1}
#'   \item{V2}{Simulated covariate V2, Binary variable with 0, 1. True time-dependent function is b(t)=sin(3*pi*t/4)}
#'   \item{V3}{Simulated covariate V2, Binary variable with 0, 1. True time-dependent function is b(t)=-1}
#'   \item{V4}{Simulated covariate V2, Binary variable with 0, 1. True time-dependent function is b(t)=(t/3)**2*exp(t/2)}
#'   \item{V5}{Simulated covariate V2, Binary variable with 0, 1. True time-dependent function is b(t)=exp(-1.5*t)}
#'   \item{event}{Simulated event variable, Binary varible with 0, 1 }
#'   \item{time}{Simulated time variable, Continous variable with non-negative value}
#'   \item{facility}{Simulated stratification variable. Have 5 levels}
#'   }
#' }
"sim_data_p5_f5"