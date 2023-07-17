# surtvep

`surtvep` is an R package for fitting penalized Newton's method 
for the time-varying effects model using mAIC, TIC and GIC as information criteria, 
in particular we span the parameters using B-spline basis functions. Utilities for carrying 
out post-estimation visualization, summarization, and inference are also provided.

# Introduction

Large-scale time-to-event data derived from national disease registries arise rapidly in medical studies. Detecting and accounting for time-varying effects is particularly important, as time-varying effects have already been reported in the clinical literature. However, there are currently no formal R packages for estimating the time-varying effects without pre-assuming the time-dependent function. Inaccurate pre-assumptions can greatly influence the estimation, leading to unreliable results. To address this issue, we developed a time-varying model using spline terms with penalization that does not require pre-assumption of the true time-dependent function, and implemented it in R.

Our package offers several benefits over traditional methods. Firstly, traditional methods for modeling time-varying survival models often rely on expanding the original data into a repeated measurement format. However, even with moderate sample sizes, this leads to a large and computationally burdensome working dataset. Our package addresses this issue by proposing a computationally efficient Kronecker product-based proximal algorithm, which allows for the evaluation of time-varying effects in large-scale studies. Additionally, our package allows for parallel computing and can handle moderate to large sample sizes more efficiently than current methods.


In our statistical software tutorial, we address a common issue encountered when analyzing data with binary covariates with near-zero variation. For example, in the SEER prostate cancer data, only 0.6% of the 716,553 patients had their tumors regional to the lymph nodes. In such cases, the associated observed information matrix of a Newton-type method may have a minimum eigenvalue close to zero and a large condition number. Inverting this nearly singular matrix can lead to numerical instability and the corresponding Newton updates may be confined within a small neighborhood of the initial value, resulting in estimates that are far from the optimal solutions. To address this problem, our proposed Proximal-Newtown method utilizes a modified Hessian matrix, which allows for accurate estimation in these scenarios.

## Models:
<table>
    <tr>
        <th>Method</th>
        <th>Description</th>
        <th>Example</th>
    </tr>
    <tr>
        <td>Newton</td>
        <td>
        Newton's method and Proximal Newton's method <a href="#references">[1]</a>.
        </td>
        <td><a href="https://um-kevinhe.github.io/surtvep/articles/surtvep.html#model-fitting">tutorial</a></td>
    </tr>
      <tr>
        <td>Newton's method with penalization</td>
        <td>
        Newton's method and Proximal Newton combined with P-spline or Smoothing-spline <a href="#references">[2]</a>.
        </td>
        <td><a href="https://um-kevinhe.github.io/surtvep/articles/surtvep.html#model-fitting">tutorial</a></td>
    </tr>
</table>

## Penalzation Coefficient Selection Methods:
<table>
    <tr>
        <th>Method</th>
        <th>Description</th>
        <th>Example</th>
    </tr>
    <tr>
        <td>mAIC</td>
        <td>
        A modified Akaki Information Criterion  <a href="#references">[1]</a>.
        </td>
        <td><a href="https://um-kevinhe.github.io/surtvep/articles/surtvep.html#model-fitting">tutorial</a></td>
    </tr>
    <tr>
        <td>TIC</td>
        <td>
        Takuchi Information Criterion  <a href="#references">[1]</a>.
        </td>
        <td><a href="https://um-kevinhe.github.io/surtvep/articles/surtvep.html#model-fitting">tutorial</a></td>
    </tr>
    <tr>
        <td>GIC</td>
        <td>
        Takuchi Information Criterion  <a href="#references">[1]</a>.
        </td>
        <td><a href="https://um-kevinhe.github.io/surtvep/articles/surtvep.html#model-fitting">tutorial</a></td>
    </tr>
    <tr>
        <td>Cross Validation</td>
        <td>
        Use cross validation to select the penalization coefficient  <a href="#references">[1]</a>.
        </td>
        <td><a href="https://um-kevinhe.github.io/surtvep/articles/surtvep.html#model-fitting">tutorial</a></td>
    </tr>
</table>

## Usage:

Here, we are using the Simulation study included in our packages as an example

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

<!---<a href="https://drive.google.com/uc?export=view&id=1ET7KIwGN6FVHtjduSNGYpIUf-ydkimIe"><img src="https://drive.google.com/uc?export=view&id=1ET7KIwGN6FVHtjduSNGYpIUf-ydkimIe" style="width: 650px; max-width: 100%; height: auto" title="Click to enlarge picture"/></a>-->

# Datasets

The SUPPORT dataset is available in the "surtvep" package. The following code will load the dataset in the form of a dataframe

    data("support")

## Simulated Datasets:

<table>
    <tr>
        <th>Dataset</th>
        <th>Size</th>
        <th>Dataset</th>
        <th>Data source</th>
    </tr>
    <tr>
        <td>simulN2kOP2</td>
        <td>25,000</td>
        <td>
        Dataset from simulation study in <a href="#references">[2]</a>.
        This is a simulation study with continuous time predictors
        </td>
        <td><a href="https://github.com/havakv/pycox/blob/master/pycox/simulations/relative_risk.py">simulN2kOP2Continuous</a>
    </tr>
    <tr>
        <td>simulN2kOP2Binary</td>
        <td>100,000</td>
        <td>
        Dataset from simulation study in <a href="#references">[12]</a>.
        This is a discrete time dataset with 1000 possible event-times.
        </td>
        <td><a href="https://github.com/havakv/pycox/tree/master/pycox/simulations/discrete_logit_hazard.py">simulN2kOP2Binary</a>
    </tr>
</table>


## Real Datasets:
<table>
    <tr>
        <th>Dataset</th>
        <th>Size</th>
        <th>Dataset</th>
        <th>Data source</th>
    </tr>
    <tr>
        <td>SUPPORT</td>
        <td></td>
        <td>
        The support dataset is a random sample of 1000 patients from Phases I & II of SUPPORT (Study to Understand Prognoses Preferences Outcomes and Risks of Treatment). This dataset is very good for learning how to fit highly nonlinear predictor effects. See 
        <a href="https://um-kevinhe.github.io/surtvep/articles/surtvep.html#model-fitting">tutorial</a></td> for preprocessing.
        </td>
        <td><a href="https://biostat.app.vumc.org/wiki/Main/SupportDesc">source</a>
    </tr>
</table>


# Installation

**Note:** *This package is still in its early stages of development, so please don't hesitate to report any problems you may experience.* 

The package only works for R 4.1.0+.

You can install 'surtvep' via:

    #Install the package, need to install the devtools packages:
    install.packages("devtools")
    require("remotes")
    remotes::install_github("UM-KevinHe/surtvep", ref = "openmp")

We recommand to start with <a href="https://um-kevinhe.github.io/surtvep/articles/surtvep.html#quick-start" target="_blank">tutorial</a>, as it provides an overview of the package's usage, including preprocessing, model training, selection of penalization parameters, and post-estimation procedures.


## Detailed tutorial

For detailed tutorial and model paramter explaination, please go to <a href="https://um-kevinhe.github.io/surtvep/index.html" target="_blank">here</a>

## Getting Help:

If you encounter any problems or bugs, please contact us at: [lfluo\@umich.edu](mailto:lfluo@umich.edu){.email}, [xuetao\@umich.edu](mailto:xuetao@umich.edu){.email}



# References

  \[1\] Wenbo Wu, Jeremy M G Taylor, Andrew F Brouwer, Lingfeng Luo, Jian Kang, Hui Jiang and Kevin He. Scalable proximal Methods for cause-specific hazard modeling with time-varying coefficients. *Lifetime Data Analysis*, 28(2):194-218, 2022. \[[paper](https://pubmed.ncbi.nlm.nih.gov/35092553/)\]

=======
# surtvep

`surtvep` is an R package for fitting Cox non-proportional hazards models with time-varying coefficients. Both unpenalized procedures (Newton and proximal Newton) and penalized procedures (P-splines and smoothing splines) are included using B-spline basis functions for estimating time-varying coefficients.
For penalized procedures, cross-validations, mAIC, TIC or GIC are implemented to select tuning parameters. Utilities for carrying out post-estimation visualization, summarization, point-wise confidence interval and hypothesis testing are also provided.

# Introduction

Large-scale time-to-event data derived from national disease registries arise rapidly in medical studies. Detecting and accounting for time-varying effects is particularly important, as time-varying effects have already been reported in the clinical literature. However, there are currently no formal R packages for estimating the time-varying effects without pre-assuming the time-dependent function. Inaccurate pre-assumptions can greatly influence the estimation, leading to unreliable results. To address this issue, we developed a time-varying model using spline terms with penalization that does not require pre-assumption of the true time-dependent function, and implemented it in R.

Our package offers several benefits over traditional methods. Firstly, traditional methods for modeling time-varying survival models often rely on expanding the original data into a repeated measurement format. However, even with moderate sample sizes, this leads to a large and computationally burdensome working dataset. Our package addresses this issue by proposing a computationally efficient Kronecker product-based proximal algorithm, which allows for the evaluation of time-varying effects in large-scale studies. Additionally, our package allows for parallel computing and can handle moderate to large sample sizes more efficiently than current methods.


In our statistical software tutorial, we address a common issue encountered when analyzing data with binary covariates with near-zero variation. For example, in the SEER prostate cancer data, only 0.6% of the 716,553 patients had their tumors regional to the lymph nodes. In such cases, the associated observed information matrix of a Newton-type method may have a minimum eigenvalue close to zero and a large condition number. Inverting this nearly singular matrix can lead to numerical instability and the corresponding Newton updates may be confined within a small neighborhood of the initial value, resulting in estimates that are far from the optimal solutions. To address this problem, our proposed Proximal-Newtown method utilizes a modified Hessian matrix, which allows for accurate estimation in these scenarios.

## Models:
<table>
    <tr>
        <th>Method</th>
        <th>Description</th>
        <th>Example</th>
    </tr>
    <tr>
        <td>Newton</td>
        <td>
        Newton's method and Proximal Newton's method <a href="#references">[1]</a>.
        </td>
        <td><a href="https://um-kevinhe.github.io/surtvep/articles/surtvep.html#model-fitting">tutorial</a></td>
    </tr>
      <tr>
        <td>Newton's method with penalization</td>
        <td>
        Newton's method and Proximal Newton combined with P-spline or Smoothing-spline <a href="#references">[2]</a>.
        </td>
        <td><a href="https://um-kevinhe.github.io/surtvep/articles/surtvep.html#model-fitting">tutorial</a></td>
    </tr>
</table>

## Penalzation Coefficient Selection Methods:
<table>
    <tr>
        <th>Method</th>
        <th>Description</th>
        <th>Example</th>
    </tr>
    <tr>
        <td>mAIC</td>
        <td>
        A modified Akaki Information Criterion  <a href="#references">[1]</a>.
        </td>
        <td><a href="https://um-kevinhe.github.io/surtvep/articles/surtvep.html#model-fitting">tutorial</a></td>
    </tr>
    <tr>
        <td>TIC</td>
        <td>
        Takuchi Information Criterion  <a href="#references">[1]</a>.
        </td>
        <td><a href="https://um-kevinhe.github.io/surtvep/articles/surtvep.html#model-fitting">tutorial</a></td>
    </tr>
    <tr>
        <td>GIC</td>
        <td>
        Takuchi Information Criterion  <a href="#references">[1]</a>.
        </td>
        <td><a href="https://um-kevinhe.github.io/surtvep/articles/surtvep.html#model-fitting">tutorial</a></td>
    </tr>
    <tr>
        <td>Cross Validation</td>
        <td>
        Use cross validation to select the penalization coefficient  <a href="#references">[1]</a>.
        </td>
        <td><a href="https://um-kevinhe.github.io/surtvep/articles/surtvep.html#model-fitting">tutorial</a></td>
    </tr>
</table>

## Usage:

Here, we are using the Simulation study included in our packages as an example

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

<!---<a href="https://drive.google.com/uc?export=view&id=1ET7KIwGN6FVHtjduSNGYpIUf-ydkimIe"><img src="https://drive.google.com/uc?export=view&id=1ET7KIwGN6FVHtjduSNGYpIUf-ydkimIe" style="width: 650px; max-width: 100%; height: auto" title="Click to enlarge picture"/></a>-->

# Datasets

The SUPPORT dataset is available in the "surtvep" package. The following code will load the dataset in the form of a dataframe

    data("support")

## Simulated Datasets:

<table>
    <tr>
        <th>Dataset</th>
        <th>Size</th>
        <th>Dataset</th>
        <!-- <th>Data source</th> -->
    </tr>
    <tr>
        <td>ExampleData</td>
        <td>4,000</td>
        <td>
        A simulated data set containing 2 continuous variables.
        </td>
        <!-- <td><a href="https://github.com/havakv/pycox/blob/master/pycox/simulations/relative_risk.py">simulN2kOP2Continuous</a> -->
    </tr>
    <tr>
        <td>ExampleDataBinary</td>
        <td>2,000</td>
        <td>
        A simulated data set containing 2 binary variables.
        </td>
        <!-- <td><a href="https://github.com/havakv/pycox/tree/master/pycox/simulations/discrete_logit_hazard.py">simulN2kOP2Binary</a> -->
    </tr>
        <tr>
        <td>StrataExample</td>
        <td>2,000</td>
        <td>
        A simulated data set containing 2 binary variables. Subjects in different strata have 
        </td>
        <!-- <td><a href="https://github.com/havakv/pycox/tree/master/pycox/simulations/discrete_logit_hazard.py">simulN2kOP2Binary</a> -->
    </tr>
</table>


## Real Datasets:
<table>
    <tr>
        <th>Dataset</th>
        <th>Size</th>
        <th>Dataset</th>
        <th>Data source</th>
    </tr>
    <tr>
        <td>SUPPORT</td>
        <td></td>
        <td>
        The support dataset is a random sample of 1000 patients from Phases I & II of SUPPORT (Study to Understand Prognoses Preferences Outcomes and Risks of Treatment). This dataset is very good for learning how to fit highly nonlinear predictor effects. See 
        <a href="https://um-kevinhe.github.io/surtvep/articles/surtvep.html#model-fitting">tutorial</a></td> for preprocessing.
        </td>
        <td><a href="https://biostat.app.vumc.org/wiki/Main/SupportDesc">source</a>
    </tr>
</table>


# Installation

**Note:** *This package is still in its early stages of development, so please don't hesitate to report any problems you may experience.* 

The package only works for R 4.1.0+.

You can install 'surtvep' via:

    #Install the package, need to install the devtools packages:
    install.packages("devtools")
    require("remotes")
    remotes::install_github("UM-KevinHe/surtvep", ref = "openmp")

We recommand to start with <a href="https://um-kevinhe.github.io/surtvep/articles/surtvep.html#quick-start" target="_blank">tutorial</a>, as it provides an overview of the package's usage, including preprocessing, model training, selection of penalization parameters, and post-estimation procedures.


## Detailed tutorial

For detailed tutorial and model paramter explaination, please go to <a href="https://um-kevinhe.github.io/surtvep/index.html" target="_blank">here</a>

## Getting Help:

If you encounter any problems or bugs, please contact us at: [lfluo\@umich.edu](mailto:lfluo@umich.edu){.email}, [xuetao\@umich.edu](mailto:xuetao@umich.edu){.email}



# References

  \[1\] Wenbo Wu, Jeremy M G Taylor, Andrew F Brouwer, Lingfeng Luo, Jian Kang, Hui Jiang and Kevin He. Scalable proximal Methods for cause-specific hazard modeling with time-varying coefficients. *Lifetime Data Analysis*, 28(2):194-218, 2022. \[[paper](https://pubmed.ncbi.nlm.nih.gov/35092553/)\]

>>>>>>> 9a2f88b98bda0338ea9c6ecf6548bc57881c8a1c
