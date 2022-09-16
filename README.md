# surtvep

"surtvep" is an R package for time-varying effects analysis in large-scale time-to-event data, including Newton's method, Proximal Newton's method and their combination with P-spline and smoothing-spline. mAIC, TIC and GIC are provided to choose the penalization coefficient.

## Get Started

Large-scale time-to-event data derived from national disease registries arise rapidly in medical studies. Detecting and accounting for time-varying effects are particularly important, as time-varying effects have already been reported in the clinical literature. However, there are no formal R packages for estimating the time-varying effect without pre-assuming the time-dependent function. However, in the real dataset, if we get the pre-assuming incorrect, the estimation could be largely influenced by the assumption. Thus, we decided to develop a time-varying model using spline terms with penalization which don't need pre-assumption for the true time-dependent function and implemented it in R.

Following are some benefits of our packages:

To begin with, traditional methods of modeling time-varying survival models typically rely on expanding the original data in a repeated measurement format, which, even with moderate sample size, usually leads to an intractably large working dataset. Consequently, the computational burden increases dramatically as the sample size grows, precluding the evaluation of time-varying effects in large-scale studies, thus, we propose a computationally efficient Kronecker product-based proximal algorithm, which enables us to extend existing methods of estimating time-varying effects for time-to-event data to a large-scale context. Detailed information about our method could be found here.

Also, by allowing parallel computing, our packages could handle the moderate and large sample size quite well compared to current methods which we would compare at the end of this page.

Specifically, when the data under analysis include a number of binary covariates with near-zero variation (e.g., in the SEER prostate cancer data, only 0.6% of the 716,553 patients had their tumors regional to the lymph nodes), the associated observed information matrix of a Newton-type method may have its minimum eigenvalue close to zero with a large condition number. Inverting such a nearly singular matrix is numerically unstable and the corresponding Newton updates are likely to be conned within a small neighborhood of the initial value, causing the estimates to be far from the optimal solutions. However, our proposed Proximal-Newtown method could handle this problem quite well by adding (how to describe the edited Hessian matrix)

## Installation

You can install 'surtvep' via:

    #Install the package, need to install the devtools packages:
    install.packages("devtools")
    require("remotes")
    remotes::install_github("UM-KevinHe/surtvep", ref = "openmp")

    #To install with Vignettes:
    install.packages("devtools")
    remotes::install_github("UM-KevinHe/surtvep", ref = "openmp", build_vignettes =T)

We recommand to start with <a href="https://um-kevinhe.github.io/surtvep/index.html" target="_blank">tutorial</a>, which explains the general usage of the package in terms of preprocessing, model training, penalziation parameter selection and evalutaion procedure.

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
        Newton's method without penalization <a href="#references">[1]</a>.
        </td>
        <td><a href="https://um-kevinhe.github.io/surtvep/articles/surtvep.html#model-fitting">tutorial</a></td>
    </tr>
    <tr>
        <td>Proximal Newton</td>
        <td>
        Proximal Newton's method without penalization <a href="#references">[2]</a>.
        </td>
        <td><a href="https://um-kevinhe.github.io/surtvep/articles/surtvep.html#model-fitting">tutorial</a></td>
    </tr>
      <tr>
        <td>Newton's method with penalization</td>
        <td>
        Proximal Newton's method without penalization <a href="#references">[2]</a>.
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

## Datasets

The SUPPORT dataset is available in the "surtvep" package. The following code will load the dataset in the form of a dataframe

    data("support")

## Example

We will use the SUPPORT dataset as an example to show the usage of "surtvep"

    event = support$death
    time  = support$d.time
    z = as.matrix(as.numeric(as.factor(support$sex)))
    colnames(z) = "gender"

    fit = coxtp(event = event, time = time, z = z)

To get the estimated time-varying effect of a specific coefficient, we could use the following plot function in our package:

    coxtp.plot(fit, coef = "gender")

## Detailed tutorial

For detailed tutorial and model paramter explaination, please go to <a href="https://um-kevinhe.github.io/surtvep/index.html" target="_blank">here</a>

## Getting Help:

If you encounter any problems or bugs, please contact us at: [lfluo\@umich.edu](mailto:lfluo@umich.edu){.email}, [xuetao\@umich.edu](mailto:xuetao@umich.edu){.email}
