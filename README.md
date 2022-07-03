# PolyMR

## Polynomial Mendelian randomization

<!-- badges: start -->

<!-- badges: end -->

This package implements the PolyMR method (Sulc et al. 2022) for estimating non-linear causal effects of an exposure on an outcome. In short, this Mendelian randomization-based method uses genotypes as instrumental variables (IVs) to model the effect of an exposure on an outcome, correcting for confounding effects using a control function. The outcome is modeled as a sum of two polynomial functions, one for the causal effect and one for the control function.

## Installation

You can install PolyMR using the <code>devtools</code> package:

``` r
devtools::install_github("JonSulc/PolyMR")
```

## Example

The main function provided by PolyMR is <code>polymr()</code>. For a dataset with n individuals and m independent IVs, the function can be called directly supplying n-length vectors for both exposure and outcome and an n-by-m matrix of genotypes.

``` r
library(PolyMR)

polymr(exposure, outcome, genotypes)
```

By default, the function uses a fixed, 10th degree polynomial for both exposure and control function, with no additional filtering of IVs. The settings used in the article differ from this in that they include filtering out IVs which have a stronger association with the outcome than the exposure and includes a model selection process in which non-significant exposure terms are iteratively eliminated until all remaining coefficients for exposure terms are significant.

We recommend the IV filtering process as any eliminated IVs would violate MR assumptions but recognize this is not necessarily a standard procedure in the field. The model selection process may result in a better fit with fewer terms, however post-selection inference bias will result in underestimated confidence intervals (or, equivalently, overestimated significance). Incorporating either of these procedures should be the result of careful consideration. The settings used in the article can be obtained by passing the following options in the method call:

``` r
polymr(exposure, outcome, genotypes, reverse_t_thr = 0, p_thr_drop = NULL)
```

A value of 0 for <code>reverse_t\_thr</code> indicates a simple comparison between effect sizes: if the variance that the IV explains in the outcome is greater than that explained in the exposure, the IV is discarded. An alternative would be to filter out only IVs that explain *significantly* more variance in the outcome than the exposure, which could be achieved by passing <code>reverse_t\_thr = qnorm(0.95)</code> (for nominal significance).

The <code>NULL</code> value passed here for <code>p_thr_drop</code> results in a Bonferroni-corrected threshold, i.e. 0.05 / number of exposure coefficients. A fixed numeric value can be passed instead.

### Simulated data

Data simulations can be performed using the (non-exported) <code>new_PolyMRDataSim()</code> function. See <code>?new_PolyMRDataSim</code>.

For example a setting with an exponential causal function (with coefficient 0.05) can be simulated with:

``` r
PolyMR:::new_PolyMRDataSim(sample_size = 10000, causal_function = function(x) 0.05*exp(x))
```
