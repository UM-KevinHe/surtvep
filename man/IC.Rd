% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IC.r
\name{IC}
\alias{IC}
\title{calculating Information Criteria from a \code{coxtp} object}
\usage{
IC(fit, IC.prox, ...)
}
\arguments{
\item{fit}{model get from coxtp}

\item{IC.prox}{when calculating information criteria, there might be numerical issue(second order derivative,
Hessian matrix approximate is singular). In such cases warnings will be given.
if \code{true}, we modified the diagonal of hessian matrix a bit, which may lead to bias issues.
Default is \code{FALSE}.}
}
\value{
\item{model.AIC}{an object with S3 class \code{"coxtp"} using AIC to select the tunning parameter.}
\item{model.TIC}{an object with S3 class \code{"coxtp"} using TIC to select the tunning parameter.}
\item{model.GIC}{an object with S3 class \code{"coxtp"} using GIC to select the tunning parameter.}
\item{AIC}{a sequence of AIC values for the different tuning parameters \code{lambda.spline} from \code{"coxtp"}.}
\item{TIC}{a sequence of AIC values for the different tuning parameters \code{lambda.spline} from \code{"coxtp"}.}
\item{GIC}{a sequence of AIC values for the different tuning parameters \code{lambda.spline} from \code{"coxtp"}.}
}
\description{
calculating Information Criteria from a \code{coxtp} object
}
\details{
In order to select the proper smoothing parameter, we utilize the idea of information criteria.
We provide four different information criteria to select the optimal smoothing parameter \eqn{\lambda}.
Generally, AIC, TIC and GIC selects similar parameters and the difference of resulting estimations are barely noticable.
See details in Lingfeng Luo et al. (2022)
}
\examples{
data(ExampleData)
z <- ExampleData$x
time <- ExampleData$time
event <- ExampleData$event
fit <- coxtp(event = event, z = z, time = time)
IC  <- IC(fit)


}