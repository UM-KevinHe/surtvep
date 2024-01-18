# surtvep

`surtvep` is an R package for fitting Cox non-proportional hazards models with time-varying coefficients.
Both unpenalized procedures (Newton and proximal Newton) and penalized procedures (P-splines and smoothing splines) are included using B-spline basis functions for estimating time-varying coefficients.
For penalized procedures, cross validations, mAIC, TIC or GIC are implemented to select tuning parameters.
Utilities for carrying out post-estimation visualization, summarization, point-wise confidence interval and hypothesis testing are also provided.

# Introduction

Large-scale time-to-event data derived from national disease registries arise rapidly in medical studies. Detecting and accounting for time-varying effects is particularly important, as time-varying effects have already been reported in the clinical literature. However, there are currently no formal R packages for estimating the time-varying effects without pre-assuming the time-dependent function. Inaccurate pre-assumptions can greatly influence the estimation, leading to unreliable results. To address this issue, we developed a time-varying model using spline terms with penalization that does not require pre-assumption of the true time-dependent function, and implemented it in R.

Our package offers several benefits over traditional methods. Firstly, traditional methods for modeling time-varying survival models often rely on expanding the original data into a repeated measurement format. However, even with moderate sample sizes, this leads to a large and computationally burdensome working dataset. Our package addresses this issue by proposing a computationally efficient Kronecker product-based proximal algorithm, which allows for the evaluation of time-varying effects in large-scale studies. Additionally, our package allows for parallel computing and can handle moderate to large sample sizes more efficiently than current methods.


In our statistical software tutorial, we address a common issue encountered when analyzing data with binary covariates with near-zero variation. For example, in the SEER prostate cancer data, only 0.6% of the 716,553 patients had their tumors regional to the lymph nodes. In such cases, the associated observed information matrix of a Newton-type method may have a minimum eigenvalue close to zero and a large condition number. Inverting this nearly singular matrix can lead to numerical instability and the corresponding Newton updates may be confined within a small neighborhood of the initial value, resulting in estimates that are far from the optimal solutions. To address this problem, our proposed Proximal-Newtown method utilizes a modified Hessian matrix, which allows for accurate estimation in these scenarios.




# Installation

**Note:** *This package is still in its early stages of development, so please don't hesitate to report any problems you may experience.* 

The package only works for R 4.1.0+.

You can install 'surtvep' via CRAN or github:

    install.packages("surtvep")

    #or
    require("devtools")
    require("remotes")
    remotes::install_github("UM-KevinHe/surtvep")

We recommand to start with <a href="https://um-kevinhe.github.io/surtvep/articles/surtvep.html#quick-start" target="_blank">tutorial</a>, as it provides an overview of the package's usage, including preprocessing, model training, selection of penalization parameters, and post-estimation procedures.


## Detailed tutorial

For detailed tutorial and model paramter explaination, please go to <a href="https://um-kevinhe.github.io/surtvep/index.html" target="_blank">here</a>.

## Getting Help:

If you encounter any problems or bugs, please contact us at: [lfluo\@umich.edu](mailto:lfluo@umich.edu){.email}, [kevinhe\@umich.edu](mailto:kevinhe@umich.edu){.email}, [Wenbo.Wu\@nyulangone.org](mailto:Wenbo.Wu@nyulangone.org){.email}

## Contributing

We welcome contributions to the **surtvep** package. Please see our [CONTRIBUTING.md](https://github.com/UM-KevinHe/surtvep/blob/main/.github/CONTRIBUTING.md) file for detailed guidelines of how to contribute.

# References

  \[1\] Gray, R. J. (1992). Flexible methods for analyzing survival data using splines, with applications to breast
cancer prognosis. *Journal of the American Statistical Association*, 87(420), 942–951. https://doi.org/10.2307/2290630

  \[2\] Gray, R. J. (1994). Spline-based tests in survival analysis. *Biometrics*, 50(3), 640–652. https://doi.org/10.2307/2532779

  \[3\] He, K., Zhu, J., Kang, J., & Li, Y. (2022). Stratified Cox models with time-varying effects for national
kidney transplant patients: A new blockwise steepest ascent method. *Biometrics*, 78(3), 1221–1232.
https://doi.org/10.1111/biom.13473

  \[4\] Luo, L., He, K., Wu, W., & Taylor, J. M. (2023). Using information criteria to select smoothing parameters when analyzing survival data with time-varying coefficient hazard models. *Statistical Methods in Medical Research*, in press.
https://doi.org/10.1177/09622802231181471

  \[5\] Wu, W., Taylor, J. M., Brouwer, A. F., Luo, L., Kang, J., Jiang, H., & He, K. (2022). Scalable proximal methods for cause-specific hazard modeling with time-varying coefficients. *Lifetime Data Analysis*, 28 (2), 194–218. https://doi.org/10.1007/s10985-021-09544-2
  
  
