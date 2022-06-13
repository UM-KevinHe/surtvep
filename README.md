
# BWslection
## Overview
BWselection is a function to do the backward selection in R based on AIC and F-test P-values. For the linear regression, the function use `lm` to call linear regression. For Logistic regression, the function use `glm` with `family="Binomial"` to call logistic regression. Currently, the function support the linear regression and the logistic regression. Just call the function by `BWselection()` with the arguments. For more resource of how to use the package, please refer to the help page and vignettes.

## Installation

```{r }
#Install the package, need to install the devtools packages:
install.packages("devtools")
devtools::install_github("xuetao666/BWselection")

#To install with Vignettes:
install.packages("devtools")
devtools::install_github("xuetao666/BWselection",build_vignettes =T)

```
## Usage:

Here, we are using the NHANES data as example:

```{r }
library(BWselection)

#Load NHANES DATA
library(NHANES)
data("NHANES")
data=NHANES

#Created variable lists for the model selection
qvarlist<-c("BPSysAve","SleepHrsNight","TotChol")
lvarlist<-c("Age","BPDiaAve","Weight","Height")
spvarlist<-c("Pulse","DirectChol")
spclist<-c(70,1.2)
catvarlist<-c("Depressed","Marijuana","Gender")

BWselection(data=NHANES,qvarlist=qvarlist,lvarlist=lvarlist,spvarlist=spvarlist,
                     spclist=spclist,catvarlist=catvarlist,outcome="BMI",type="lm",sig=0.05,complete_case=TRUE)

```
## Improvement compared to step function in R

* 1. In the long-term list of variables needed to be selected from, this version only need to input the variable list instead of the formalized formula. For example, in the step version, one need to fit the model first, then do the selection, in the model fit step, formulas like "A~B+C+D+E+..." needed to be writened down carefully, including the split terms. If one wants to change the orginal model, it is hard to modifiy and could higher chances to make mistakes. However, in our model, you only need to specify the variable type and list in the function, the function will generate the formula for you.
* 2. This version handles quandratic terms better. The "step" function might remove the linear term of quandratic form and leave the quandratic term in which is not correct

## Getting Help:
If you encounter any problems or bugs, please contact me at xuetao@umich.edu

