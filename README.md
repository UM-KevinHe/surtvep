
# surtvep
Cox Non-PH model with penalization

## Overview
Large-scale time-to-event data derived from national disease registries arise rapidly in medical studies. Detecting and accounting for time-varying effects are particularly important, as time-varying effects have already been reported in the clinical literature. However, there are no formal R packages for estimating the time-varying effect without pre-assuming the time-dependent function. However, in the real dataset, if we get the pre-assuming incorrect, the estimation could be largely influenced by the assumption. Thus, we decided to develop a time-varying model using spline terms with penalization which don't need pre-assumption for the true time-dependent function and implemented it in R.

Following are some benefits of our packages: 

To begin with, traditional methods of modeling time-varying survival models typically rely on expanding the original data in a repeated measurement format, which, even with moderate sample size, usually leads to an intractably large working dataset. Consequently, the computational burden increases dramatically as the sample size grows, precluding the evaluation of time-varying effects in large-scale studies, thus, we propose a computationally efficient Kronecker product-based proximal algorithm, which enables us to extend existing methods of estimating time-varying effects for time-to-event data to a large-scale context. Detailed information about our method could be found here. 

Also, by allowing parallel computing, our packages could handle the moderate and large sample size quite well compared to current methods which we would compare at the end of this page. 

Specifically, when the data under analysis include a number of binary covariates with near-zero variation (e.g., in the SEER prostate cancer data, only 0.6% of the 716,553 patients had their tumors regional to the lymph nodes), the associated observed information matrix of a Newton-type method may have its minimum eigenvalue close to zero with a large condition number. Inverting such a nearly singular matrix is numerically unstable and the corresponding Newton updates are likely to be conned within a small neighborhood of the initial value, causing the estimates to be far from the optimal solutions. However, our proposed Proximal-Newtown method could handle this problem quite well by adding (how to describe the edited Hessian matrix)


## Installation

```{r }
#Install the package, need to install the devtools packages:
install.packages("devtools")
devtools::install_github("UM-KevinHe/surtvep")

#To install with Vignettes:
install.packages("devtools")
devtools::install_github("UM-KevinHe/surtvep",build_vignettes =T)

```
## Usage:

Here, we are using the Simulation study included in our packages as an example

```{r }
library(surtvep)

#Load Simulation study
sim_data=sim_data
#Clean and create label and covariate matrix for the package:
event=sim_data[,"event"]
time=sim_data[,"time"]
data=sim_data[,!colnames(sim_data) %in% c("event","time")]

#Fit the model(Time varying model without penalty)

fit <- coxtp(event = event, z = data, time = time)
coxtp.plot(fit,coef="V1")

```
<a href="https://drive.google.com/uc?export=view&id=1ET7KIwGN6FVHtjduSNGYpIUf-ydkimIe"><img src="https://drive.google.com/uc?export=view&id=1ET7KIwGN6FVHtjduSNGYpIUf-ydkimIe" style="width: 650px; max-width: 100%; height: auto" title="Click to enlarge picture"/></a>

  
## Penalized Newton's method
In order to add penalty to the penalized Newton's method, three additional parameters can be added.

Now, set a sequence of smoothing paramemter to choose from. The default range is set as
```{r example.fit, eval=FALSE}
#specify smoothing parameter lambda:
lambda_spline = c(1:10)

#if not specified, the default range is set as follows, where n is the sample size
lambda_spline = seq(0, n^(1/3), length.out = 10)
```

Then set the relevant parameters and fit a model to the prepared data:
```{r example.fit, eval=FALSE}
fit.spline <- surtvep(Surv(time, event)~., data = simulData, 
               lambda = c(0:100),
               spline = "Smooth-spline",
               IC = "all")
}
```



## Detailed tutorial

  
For detailed tutorial and model paramter explaination, please go to <a href="https://um-kevinhe.github.io/surtvep/index.html" target="_blank">here</a>

## Getting Help:
If you encounter any problems or bugs, please contact me at:  xuetao@umich.edu

