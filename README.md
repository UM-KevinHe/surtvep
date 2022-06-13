
# surtvep
Cox Non-PH model with penalization

## Overview
Large-scale time-to-event data derived from national disease registries arise rapidly in medical studies. Detecting and accounting for time-varying effects are particularly important, as time-varying effects have already been reported in the clinical literature.  However, traditional methods of modeling time-varying survival models typically rely on expanding the original data in a repeated measurement format, which, even with moderate sample size, usually leads to an intractably large working dataset. Consequently, the computational burden increases dramatically as the sample size grows, precluding the evaluation of time-varying effects in large-scale studies.  

Thus, propose a computationally efficient Kronecker product-based proximal algorithm, which enables us to extend existing methods of estimating time-varying effects for time-to-event data to a large-scale context. Detailed information about our method could be found under (Insert Link) Also, by allowing paralle computing, our packages could handle the moderate and large sample size quite well compared to current packages and methods which we would compare at the end of this file.


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
## Improvement compared to previous Time-varying packages:



## Getting Help:
If you encounter any problems or bugs, please contact me at xuetao@umich.edu

