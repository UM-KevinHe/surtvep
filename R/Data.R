#' Example data for surtvep with 4000 observations of 2 variables:
#' 
#' A simulated data set containing two 2 variables. 
#' @name ExampleData
#' @docType data
#' @usage data(ExampleData)
#'
#' @format A List containing the following elements:
#' \describe{
#'   \item{x}{Simulated covariate V1 and v2, Continuous variable. True time-dependent function is b(t)=1 and true time-dependent function is b(t)=sin(3*pi*t/4)}
#'   \item{event}{Simulated event variable, Binary varible with 0, 1 }
#'   \item{time}{Simulated time variable, Continous variable with non-negative value}
#' }
"ExampleData"

#' ExampleDataBinary
#' Example data for surtvep with 2000 observations of 2 binary variables:
#' 
#' A simulated data set containing two 2 variables. 
#' @name ExampleDataBinary
#' @docType data
#' @usage data(ExampleDataBinary)
#'
#' @format A List containing the following elements:
#' \describe{
#'   \item{V1}{Simulated covariate V1, Binary variable. True time-dependent function is b(t)=1}
#'   \item{V2}{Simulated covariate V2, Binary variable. True time-dependent function is b(t)=exp(-1.5*t)}
#'   \item{event}{Simulated event variable, Binary varible with 0, 1 }
#'   \item{time}{Simulated time variable, Continous variable with non-negative value}
#' }
"ExampleDataBinary"