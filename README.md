
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`InterXShift` <img src="man/figures/InterXshift_logo.png" style="float:right; height:200px;">

<!-- badges: start -->

[![R-CMD-check](https://github.com/blind-contours/InterXShift/workflows/R-CMD-check/badge.svg)](https://github.com/blind-contours/InterXShift/actions)
[![Coverage
Status](https://img.shields.io/codecov/c/github/blind-contours/InterXShift/master.svg)](https://codecov.io/github/blind-contours/InterXShift?branch=master)
[![CRAN](https://www.r-pkg.org/badges/version/InterXShift)](https://www.r-pkg.org/pkg/InterXShift)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/InterXShift)](https://CRAN.R-project.org/package=InterXShift)
[![CRAN total
downloads](http://cranlogs.r-pkg.org/badges/grand-total/InterXShift)](https://CRAN.R-project.org/package=InterXShift)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![MIT
license](https://img.shields.io/badge/license-MIT-brightgreen.svg)](https://opensource.org/licenses/MIT)
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4070042.svg)](https://doi.org/10.5281/zenodo.4070042) -->
<!-- [![DOI](https://joss.theoj.org/papers/10.21105/joss.02447/status.svg)](https://doi.org/10.21105/joss.02447) -->
<!-- badges: end -->

> Interaction Identification and Estimation using Data-Adaptive
> Stochastic Interventions **Authors:** [David
> McCoy](https://davidmccoy.org)

------------------------------------------------------------------------

## What’s `InterXshift`?

The `InterXshift` R package offers an approach which identifies and
estimates the impact interactions in a mixed exposure on an outcome. We
define interaction as the counterfactual mean of the outcome under
stochastic interventions of two exposures compared to the additive
counterfactual mean of the two expowures intervened on indepdentently.
These interventions or exposure changes depend on naturally observed
values, as described in past literature (Dı́az and van der Laan 2012;
Haneuse and Rotnitzky 2013).

`InterXshift` builds on work described in (McCoy et al. 2023). However
instead of identifying interactions through an semi-parametric
definition of an F-statistics and then estimating our interaction target
parameter using CV-TMLE pooled across exposure sets, we provide a more
streamlined approach. In this package, we identify interactions through
g-computation first - evaluating the expected outcome under joint shift
compared to the sum of individual shifts using Super Learner. We then
rank these estimates as the highest sets of synergy and antagonism. We
then use CV-TMLE and pool within the ranks.

The package ensures robustness by employing a k-fold cross-validation
framework. This framework helps in estimating a data-adaptive parameter,
which is the stochastic shift target parameters for the exposure sets
identified as having synergy or antagonism. The process begins by
partitioning the data into parameter-generating and estimation samples.
In the parameter-generating sample, we identify our ranks of
antagonistic and synergistic exposure sets through a machine learning
g-computation framework. In the estimation sample we then estimate our
interaction target parameter using the doubly robust estimator TMLE to
ensure asymptotic efficiency which allows us to construct confidence
intervals for our estimates (unlike the g-comp method).

By using InterXshift, users get access to a tool that offers both k-fold
specific and aggregated results for the top synergistic and antagonistic
relationships, ensuring that researchers can glean the most information
from their data. For a more in-depth exploration, there’s an
accompanying vignette.

To utilize the package, users need to provide vectors for exposures,
covariates, and outcomes. They also specify the respective $\delta$ for
each exposure (indicating the degree of shift) and if this delta should
be adaptive in response to positivity violations. The `top_n` parameter
defines the top number of synergistic, antagonistic, positive and
negative ranked impacts to estiamte. A detailed guide is provided in the
vignette. With these inputs, `InterXshift` processes the data and
delivers tables showcasing fold-specific results and aggregated
outcomes, allowing users to glean insights effectively.

`InterXshift` also incorporates features from the `sl3` package (Coyle,
Hejazi, Malenica, et al. 2022), facilitating ensemble machine learning
in the estimation process. If the user does not specify any stack
parameters, `InterXshift` will automatically create an ensemble of
machine learning algorithms that strike a balance between flexibility
and computational efficiency.

------------------------------------------------------------------------

## Installation

*Note:* Because the `InterXshift` package (currently) depends on `sl3`
that allows ensemble machine learning to be used for nuisance parameter
estimation and `sl3` is not on CRAN the `InterXshift` package is not
available on CRAN and must be downloaded here.

There are many depedencies for `InterXshift` so it’s easier to break up
installation of the various packages to ensure proper installation.

First install the basis estimators used in the data-adaptive variable
discovery of the exposure and covariate space:

`InterXshift` uses the `sl3` package to build ensemble machine learners
for each nuisance parameter.

``` r
remotes::install_github("tlverse/sl3@devel")
```

Make sure `sl3` installs correctly then install `InterXshift`

``` r
remotes::install_github("blind-contours/InterXshift@main")
```

------------------------------------------------------------------------

## Example

To illustrate how `InterXshift` may be used to ascertain the effect of a
mixed exposure, consider the following example:

``` r
library(InterXshift)
library(devtools)
#> Loading required package: usethis
library(kableExtra)
library(sl3)

seed <- 429153
set.seed(seed)
```

We will directly use synthetic data from the NIEHS used to test new
mixture methods. This data has built in strong positive and negative
marginal effects and certain interactions. Found here:
<https://github.com/niehs-prime/2015-NIEHS-MIxtures-Workshop>

``` r
data("NIEHS_data_1", package = "InterXshift")
```

``` r
NIEHS_data_1$W <- rnorm(nrow(NIEHS_data_1), mean = 0, sd = 0.1)
w <- NIEHS_data_1[, c("W", "Z")]
a <- NIEHS_data_1[, c("X1", "X2", "X3", "X4", "X5", "X6", "X7")]
y <- NIEHS_data_1$Y

deltas <- list(
  "X1" = 1, "X2" = 1, "X3" = 1,
  "X4" = 1, "X5" = 1, "X6" = 1, "X7" = 1
)
head(NIEHS_data_1) %>%
  kbl(caption = "NIEHS Data") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
NIEHS Data
</caption>
<thead>
<tr>
<th style="text-align:right;">
obs
</th>
<th style="text-align:right;">
Y
</th>
<th style="text-align:right;">
X1
</th>
<th style="text-align:right;">
X2
</th>
<th style="text-align:right;">
X3
</th>
<th style="text-align:right;">
X4
</th>
<th style="text-align:right;">
X5
</th>
<th style="text-align:right;">
X6
</th>
<th style="text-align:right;">
X7
</th>
<th style="text-align:right;">
Z
</th>
<th style="text-align:right;">
W
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
7.534686
</td>
<td style="text-align:right;">
0.4157066
</td>
<td style="text-align:right;">
0.5308077
</td>
<td style="text-align:right;">
0.2223965
</td>
<td style="text-align:right;">
1.1592634
</td>
<td style="text-align:right;">
2.4577556
</td>
<td style="text-align:right;">
0.9438601
</td>
<td style="text-align:right;">
1.8714406
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.1335790
</td>
</tr>
<tr>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
19.611934
</td>
<td style="text-align:right;">
0.5293572
</td>
<td style="text-align:right;">
0.9339570
</td>
<td style="text-align:right;">
1.1210595
</td>
<td style="text-align:right;">
1.3350074
</td>
<td style="text-align:right;">
0.3096883
</td>
<td style="text-align:right;">
0.5190970
</td>
<td style="text-align:right;">
0.2418065
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0585291
</td>
</tr>
<tr>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
12.664050
</td>
<td style="text-align:right;">
0.4849759
</td>
<td style="text-align:right;">
0.7210988
</td>
<td style="text-align:right;">
0.4629027
</td>
<td style="text-align:right;">
1.0334138
</td>
<td style="text-align:right;">
0.9492810
</td>
<td style="text-align:right;">
0.3664090
</td>
<td style="text-align:right;">
0.3502445
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.1342057
</td>
</tr>
<tr>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
15.600288
</td>
<td style="text-align:right;">
0.8275456
</td>
<td style="text-align:right;">
1.0457137
</td>
<td style="text-align:right;">
0.9699040
</td>
<td style="text-align:right;">
0.9045099
</td>
<td style="text-align:right;">
0.9107914
</td>
<td style="text-align:right;">
0.4299847
</td>
<td style="text-align:right;">
1.0007901
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0734320
</td>
</tr>
<tr>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
18.606498
</td>
<td style="text-align:right;">
0.5190363
</td>
<td style="text-align:right;">
0.7802400
</td>
<td style="text-align:right;">
0.6142188
</td>
<td style="text-align:right;">
0.3729743
</td>
<td style="text-align:right;">
0.5038126
</td>
<td style="text-align:right;">
0.3575472
</td>
<td style="text-align:right;">
0.5906156
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
-0.0148427
</td>
</tr>
<tr>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
18.525890
</td>
<td style="text-align:right;">
0.4009491
</td>
<td style="text-align:right;">
0.8639886
</td>
<td style="text-align:right;">
0.5501847
</td>
<td style="text-align:right;">
0.9011016
</td>
<td style="text-align:right;">
1.2907615
</td>
<td style="text-align:right;">
0.7990418
</td>
<td style="text-align:right;">
1.5097039
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.1749775
</td>
</tr>
</tbody>
</table>

Based on the data key, we expect X1 to have the strongest positive
effect, X5 the strongest negative. So we would expect these to take the
top ranks for these marginal associations. For interactions we expect to
find interactions built into the data if there is enough observations in
the region. The github page gives details on the types of interactions
built into this synthetic data.

``` r

ptm <- proc.time()
sim_results <- InterXshift(
  w = w,
  a = a,
  y = y,
  delta = deltas,
  n_folds = 3,
  num_cores = 6,
  outcome_type = "continuous",
  seed = seed,
  top_n = 2
)
#> 
#> Iter: 1 fn: 188.2217  Pars:  0.18996 0.81004
#> Iter: 2 fn: 188.2217  Pars:  0.18996 0.81004
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 376.1549  Pars:  0.15230 0.84770
#> Iter: 2 fn: 376.1549  Pars:  0.15230 0.84770
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 330.7270  Pars:  0.0000001547 0.9999998440
#> Iter: 2 fn: 330.7270  Pars:  0.00000009486 0.99999990514
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 395.5675  Pars:  0.9999995911 0.0000004086
#> Iter: 2 fn: 395.5675  Pars:  0.9999997531 0.0000002469
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 329.0977  Pars:  0.00000009094 0.99999990921
#> Iter: 2 fn: 329.0977  Pars:  0.00000005678 0.99999994322
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 394.9951  Pars:  0.97308 0.02692
#> Iter: 2 fn: 394.9951  Pars:  0.97315 0.02685
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 411.6276  Pars:  0.77643 0.22357
#> Iter: 2 fn: 411.6276  Pars:  0.77645 0.22355
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 396.0738  Pars:  0.9999991806 0.0000008196
#> Iter: 2 fn: 396.0738  Pars:  0.9999995607 0.0000004393
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 359.4551  Pars:  0.001522 0.998478
#> Iter: 2 fn: 358.7028  Pars:  0.57483 0.42517
#> Iter: 3 fn: 358.7028  Pars:  0.57483 0.42517
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 271.4984  Pars:  0.05346 0.94654
#> Iter: 2 fn: 271.4984  Pars:  0.05345 0.94655
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 138.6347  Pars:  0.14848 0.85152
#> Iter: 2 fn: 138.6347  Pars:  0.14848 0.85152
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 66.5503   Pars:  0.0000001271 0.9999998736
#> Iter: 2 fn: 66.5503   Pars:  0.00000007528 0.99999992472
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: -5.0555   Pars:  0.44462 0.55538
#> Iter: 2 fn: -5.0555   Pars:  0.44462 0.55538
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 183.8985  Pars:  0.20871 0.79129
#> Iter: 2 fn: 183.8985  Pars:  0.20871 0.79129
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 383.7304  Pars:  0.15233 0.84767
#> Iter: 2 fn: 383.7304  Pars:  0.15232 0.84768
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 380.2020  Pars:  0.23886 0.76114
#> Iter: 2 fn: 380.2020  Pars:  0.23886 0.76114
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 183.6422  Pars:  0.03878 0.96122
#> Iter: 2 fn: 183.6422  Pars:  0.03878 0.96122
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 333.7236  Pars:  0.03373 0.96627
#> Iter: 2 fn: 333.7236  Pars:  0.03373 0.96627
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 341.5709  Pars:  0.21640 0.78360
#> Iter: 2 fn: 341.5709  Pars:  0.21640 0.78360
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 388.3596  Pars:  0.34865 0.65135
#> Iter: 2 fn: 388.3596  Pars:  0.34865 0.65135
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 342.6278  Pars:  0.60993 0.39007
#> Iter: 2 fn: 342.6278  Pars:  0.60876 0.39124
#> Iter: 3 fn: 342.6278  Pars:  0.60876 0.39124
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 388.2628  Pars:  0.20855 0.79145
#> Iter: 2 fn: 388.2628  Pars:  0.20855 0.79145
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 393.7808  Pars:  0.58218 0.41782
#> Iter: 2 fn: 393.7808  Pars:  0.58218 0.41782
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 390.2628  Pars:  0.19613 0.80387
#> Iter: 2 fn: 390.2628  Pars:  0.19612 0.80388
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 360.8686  Pars:  0.999997556 0.000002444
#> Iter: 2 fn: 360.8686  Pars:  0.999998373 0.000001627
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 299.4745  Pars:  0.20485 0.79515
#> Iter: 2 fn: 299.4745  Pars:  0.20486 0.79514
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 185.8446  Pars:  0.03814 0.96186
#> Iter: 2 fn: 185.8446  Pars:  0.03814 0.96186
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 86.1185   Pars:  0.26299 0.73701
#> Iter: 2 fn: 86.1185   Pars:  0.26298 0.73702
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 18.4205   Pars:  0.37690 0.62310
#> Iter: 2 fn: 18.4205   Pars:  0.37690 0.62310
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 184.0532  Pars:  0.00000009628 0.99999990364
#> Iter: 2 fn: 184.0532  Pars:  0.00000005058 0.99999994942
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 356.5325  Pars:  0.18829 0.81171
#> Iter: 2 fn: 356.5325  Pars:  0.18829 0.81171
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 353.1351  Pars:  0.24568 0.75432
#> Iter: 2 fn: 353.1351  Pars:  0.24568 0.75432
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 180.6660  Pars:  0.00000004153 0.99999995880
#> Iter: 2 fn: 180.6660  Pars:  0.00000002203 0.99999997797
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 390.6584  Pars:  0.04254 0.95746
#> Iter: 2 fn: 390.6584  Pars:  0.04253 0.95747
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 345.8211  Pars:  0.03422 0.96578
#> Iter: 2 fn: 345.8211  Pars:  0.03422 0.96578
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 400.1339  Pars:  0.999997412 0.000002589
#> Iter: 2 fn: 400.1339  Pars:  0.9999993596 0.0000006404
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 347.2235  Pars:  0.000007085 0.999992915
#> Iter: 2 fn: 347.2235  Pars:  0.000001319 0.999998681
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 401.4868  Pars:  0.999997153 0.000002846
#> Iter: 2 fn: 401.4868  Pars:  0.999998259 0.000001741
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 403.9855  Pars:  0.75152 0.24848
#> Iter: 2 fn: 403.9855  Pars:  0.75154 0.24846
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 401.3947  Pars:  0.999984 0.000016
#> Iter: 2 fn: 401.3947  Pars:  0.999995323 0.000004677
#> Iter: 3 fn: 401.3947  Pars:  0.999998239 0.000001761
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 350.1500  Pars:  0.999997608 0.000002392
#> Iter: 2 fn: 350.1500  Pars:  0.99999857 0.00000143
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 256.1493  Pars:  0.09862 0.90138
#> Iter: 2 fn: 256.1493  Pars:  0.09862 0.90138
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 176.8860  Pars:  0.0000000001619 1.0000000003752
#> Iter: 2 fn: 176.8860  Pars:  4.409e-11 1.000e+00
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 83.6283   Pars:  0.14539 0.85461
#> Iter: 2 fn: 83.6283   Pars:  0.14539 0.85461
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 6.8251    Pars:  0.008576 0.991424
#> Iter: 2 fn: 6.8251    Pars:  0.00852 0.99148
#> Iter: 3 fn: 6.8251    Pars:  0.00852 0.99148
#> solnp--> Completed in 3 iterations
#> 
#> Iter: 1 fn: 192.2133  Pars:  0.05664 0.94336
#> Iter: 2 fn: 192.2133  Pars:  0.05664 0.94336
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 396.6718  Pars:  0.18042 0.81958
#> Iter: 2 fn: 396.6718  Pars:  0.18042 0.81958
#> solnp--> Completed in 2 iterations
#> 
#> Iter: 1 fn: 387.5808  Pars:  0.08945 0.91055
#> Iter: 2 fn: 387.5808  Pars:  0.08945 0.91055
#> solnp--> Completed in 2 iterations
proc.time() - ptm
#>     user   system  elapsed 
#>   67.194    3.702 1277.583

## marginal effects
top_positive_effects <- sim_results$`Pos Shift Results by Rank`
top_negative_effects <- sim_results$`Neg Shift Results by Rank`

## interaction effects
pooled_synergy_effects <- sim_results$`Pooled Synergy Results by Rank`
pooled_antagonism_effects <- sim_results$`Pooled Antagonism Results by Rank`

k_fold_synergy_effects <- sim_results$`K Fold Synergy Results`
k_fold_antagonism_effects <- sim_results$`K Fold Antagonism Results`
```

``` r
top_positive_effects$`Rank 1` %>%
  kbl(caption = "Rank 1 Positive Stochastic Intervention Results") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Rank 1 Positive Stochastic Intervention Results
</caption>
<thead>
<tr>
<th style="text-align:left;">
Condition
</th>
<th style="text-align:right;">
Psi
</th>
<th style="text-align:right;">
Variance
</th>
<th style="text-align:right;">
SE
</th>
<th style="text-align:right;">
Lower CI
</th>
<th style="text-align:right;">
Upper CI
</th>
<th style="text-align:right;">
P-value
</th>
<th style="text-align:left;">
Fold
</th>
<th style="text-align:left;">
Type
</th>
<th style="text-align:left;">
Variables
</th>
<th style="text-align:right;">
N
</th>
<th style="text-align:right;">
Delta
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
X1
</td>
<td style="text-align:right;">
13.81500
</td>
<td style="text-align:right;">
0.4971182
</td>
<td style="text-align:right;">
0.7050661
</td>
<td style="text-align:right;">
12.4331
</td>
<td style="text-align:right;">
15.1969
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
Indiv Shift
</td>
<td style="text-align:left;">
X1
</td>
<td style="text-align:right;">
167
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
X1
</td>
<td style="text-align:right;">
12.00894
</td>
<td style="text-align:right;">
0.7999055
</td>
<td style="text-align:right;">
0.8943744
</td>
<td style="text-align:right;">
10.2560
</td>
<td style="text-align:right;">
13.7619
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
Indiv Shift
</td>
<td style="text-align:left;">
X1
</td>
<td style="text-align:right;">
167
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
X1
</td>
<td style="text-align:right;">
16.34802
</td>
<td style="text-align:right;">
0.5371240
</td>
<td style="text-align:right;">
0.7328874
</td>
<td style="text-align:right;">
14.9116
</td>
<td style="text-align:right;">
17.7844
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
3
</td>
<td style="text-align:left;">
Indiv Shift
</td>
<td style="text-align:left;">
X1
</td>
<td style="text-align:right;">
166
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Rank 1
</td>
<td style="text-align:right;">
14.98981
</td>
<td style="text-align:right;">
0.2786728
</td>
<td style="text-align:right;">
0.5278947
</td>
<td style="text-align:right;">
13.9552
</td>
<td style="text-align:right;">
16.0245
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:left;">
Indiv Shift
</td>
<td style="text-align:left;">
Rank 1
</td>
<td style="text-align:right;">
500
</td>
<td style="text-align:right;">
1
</td>
</tr>
</tbody>
</table>

Above we show the findings for the top rank positive marginal effect.
Here we consistently find X1 which is true based on what is built into
the DGP. The pooled estimate is pooling the findings for the top ranked
positive result found across the folds which is all X1 in this case.

Next we look at the top negative result:

``` r
top_negative_effects$`Rank 2` %>%
  kbl(caption = "Rank 1 Negative Stochastic Intervention Results") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Rank 1 Negative Stochastic Intervention Results
</caption>
<thead>
<tr>
<th style="text-align:left;">
Condition
</th>
<th style="text-align:right;">
Psi
</th>
<th style="text-align:right;">
Variance
</th>
<th style="text-align:right;">
SE
</th>
<th style="text-align:right;">
Lower CI
</th>
<th style="text-align:right;">
Upper CI
</th>
<th style="text-align:right;">
P-value
</th>
<th style="text-align:left;">
Fold
</th>
<th style="text-align:left;">
Type
</th>
<th style="text-align:left;">
Variables
</th>
<th style="text-align:right;">
N
</th>
<th style="text-align:right;">
Delta
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
X5
</td>
<td style="text-align:right;">
-3.735474
</td>
<td style="text-align:right;">
0.6972921
</td>
<td style="text-align:right;">
0.8350402
</td>
<td style="text-align:right;">
-5.3721
</td>
<td style="text-align:right;">
-2.0988
</td>
<td style="text-align:right;">
7.7e-06
</td>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
Indiv Shift
</td>
<td style="text-align:left;">
X5
</td>
<td style="text-align:right;">
167
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
X5
</td>
<td style="text-align:right;">
-3.708447
</td>
<td style="text-align:right;">
0.6808466
</td>
<td style="text-align:right;">
0.8251343
</td>
<td style="text-align:right;">
-5.3257
</td>
<td style="text-align:right;">
-2.0912
</td>
<td style="text-align:right;">
7.0e-06
</td>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
Indiv Shift
</td>
<td style="text-align:left;">
X5
</td>
<td style="text-align:right;">
167
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
X5
</td>
<td style="text-align:right;">
-3.337799
</td>
<td style="text-align:right;">
0.4530761
</td>
<td style="text-align:right;">
0.6731093
</td>
<td style="text-align:right;">
-4.6571
</td>
<td style="text-align:right;">
-2.0185
</td>
<td style="text-align:right;">
7.0e-07
</td>
<td style="text-align:left;">
3
</td>
<td style="text-align:left;">
Indiv Shift
</td>
<td style="text-align:left;">
X5
</td>
<td style="text-align:right;">
166
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Rank 2
</td>
<td style="text-align:right;">
-3.502256
</td>
<td style="text-align:right;">
0.2273206
</td>
<td style="text-align:right;">
0.4767815
</td>
<td style="text-align:right;">
-4.4367
</td>
<td style="text-align:right;">
-2.5678
</td>
<td style="text-align:right;">
0.0e+00
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:left;">
Indiv Shift
</td>
<td style="text-align:left;">
Rank 2
</td>
<td style="text-align:right;">
500
</td>
<td style="text-align:right;">
1
</td>
</tr>
</tbody>
</table>

Here we consistently see X5 as having the strongest negative impact
which is also true compared to the true DGP.

Next we will look at the top synergy results which is defined as the
exposures that when shifted jointly have the highest, most positive,
expected outcome difference compared to the sum of individual shifts of
the same variables.

``` r
pooled_synergy_effects$`Rank 1` %>%
  kbl(caption = "Rank 1 Synergy Stochastic Intervention Results") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Rank 1 Synergy Stochastic Intervention Results
</caption>
<thead>
<tr>
<th style="text-align:left;">
Rank
</th>
<th style="text-align:right;">
Psi
</th>
<th style="text-align:right;">
Variance
</th>
<th style="text-align:right;">
SE
</th>
<th style="text-align:right;">
Lower CI
</th>
<th style="text-align:right;">
Upper CI
</th>
<th style="text-align:right;">
P-value
</th>
<th style="text-align:left;">
Fold
</th>
<th style="text-align:right;">
N
</th>
<th style="text-align:right;">
Delta Exposure 1
</th>
<th style="text-align:right;">
Delta Exposure 2
</th>
<th style="text-align:left;">
Type
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Rank 1
</td>
<td style="text-align:right;">
-0.7658578
</td>
<td style="text-align:right;">
0.2417480
</td>
<td style="text-align:right;">
0.4916788
</td>
<td style="text-align:right;">
-1.7295
</td>
<td style="text-align:right;">
0.1978
</td>
<td style="text-align:right;">
0.2747394
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:right;">
500
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Var 1
</td>
</tr>
<tr>
<td style="text-align:left;">
Rank 1
</td>
<td style="text-align:right;">
-3.5759161
</td>
<td style="text-align:right;">
0.2393057
</td>
<td style="text-align:right;">
0.4891888
</td>
<td style="text-align:right;">
-4.5347
</td>
<td style="text-align:right;">
-2.6171
</td>
<td style="text-align:right;">
0.0000003
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:right;">
500
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Var 2
</td>
</tr>
<tr>
<td style="text-align:left;">
Rank 1
</td>
<td style="text-align:right;">
-3.9894164
</td>
<td style="text-align:right;">
0.2265181
</td>
<td style="text-align:right;">
0.4759392
</td>
<td style="text-align:right;">
-4.9222
</td>
<td style="text-align:right;">
-3.0566
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:right;">
500
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Joint
</td>
</tr>
<tr>
<td style="text-align:left;">
Rank 1
</td>
<td style="text-align:right;">
0.3523575
</td>
<td style="text-align:right;">
0.2508620
</td>
<td style="text-align:right;">
0.5008613
</td>
<td style="text-align:right;">
-0.6293
</td>
<td style="text-align:right;">
1.3340
</td>
<td style="text-align:right;">
0.6185685
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:right;">
500
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Interaction
</td>
</tr>
</tbody>
</table>

Above this table shows the pooled results for the rank 1 synergy
exposure interaction. Of course, the exposure sets in the interaction
deemed to have the highest impact, synergy, may differ between the folds
and thus this pooling may be over different exposure sets. Thus, the
first line shows the pooled estimate for a shift in the first variable,
the second line the second variable, third line the joint and fourth
line the difference between the joint and sum of the first two lines, or
the interaction effect. Therefore, in this case, we could be pooling
over different variables because of inconcistency in what is included as
rank 1 between the folds. Next we look at the k-fold specific results.

``` r
k_fold_synergy_effects$`Rank 1` %>%
  kbl(caption = "K-fold Synergy Stochastic Intervention Results") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
K-fold Synergy Stochastic Intervention Results
</caption>
<thead>
<tr>
<th style="text-align:right;">
Rank
</th>
<th style="text-align:right;">
Psi
</th>
<th style="text-align:right;">
Variance
</th>
<th style="text-align:right;">
SE
</th>
<th style="text-align:right;">
Lower CI
</th>
<th style="text-align:right;">
Upper CI
</th>
<th style="text-align:right;">
P-value
</th>
<th style="text-align:right;">
Fold
</th>
<th style="text-align:right;">
N
</th>
<th style="text-align:right;">
Delta Exposure 1
</th>
<th style="text-align:right;">
Delta Exposure 2
</th>
<th style="text-align:left;">
Type
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.8485023
</td>
<td style="text-align:right;">
0.7309842
</td>
<td style="text-align:right;">
0.8549761
</td>
<td style="text-align:right;">
-2.5242
</td>
<td style="text-align:right;">
0.8272
</td>
<td style="text-align:right;">
0.3588033
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
167
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
X4
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-3.8157757
</td>
<td style="text-align:right;">
0.7073221
</td>
<td style="text-align:right;">
0.8410244
</td>
<td style="text-align:right;">
-5.4642
</td>
<td style="text-align:right;">
-2.1674
</td>
<td style="text-align:right;">
0.0000317
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
167
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
X5
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-3.8404800
</td>
<td style="text-align:right;">
0.6566876
</td>
<td style="text-align:right;">
0.8103627
</td>
<td style="text-align:right;">
-5.4288
</td>
<td style="text-align:right;">
-2.2522
</td>
<td style="text-align:right;">
0.0000199
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
167
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
X4-X5
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.8237980
</td>
<td style="text-align:right;">
0.9055973
</td>
<td style="text-align:right;">
0.9516288
</td>
<td style="text-align:right;">
-1.0414
</td>
<td style="text-align:right;">
2.6890
</td>
<td style="text-align:right;">
0.3984039
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
167
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Interaction
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.6330329
</td>
<td style="text-align:right;">
0.6895704
</td>
<td style="text-align:right;">
0.8304037
</td>
<td style="text-align:right;">
-2.2606
</td>
<td style="text-align:right;">
0.9945
</td>
<td style="text-align:right;">
0.4872591
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
167
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
X4
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-3.7729863
</td>
<td style="text-align:right;">
0.7084949
</td>
<td style="text-align:right;">
0.8417214
</td>
<td style="text-align:right;">
-5.4227
</td>
<td style="text-align:right;">
-2.1232
</td>
<td style="text-align:right;">
0.0000391
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
167
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
X5
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-4.2714349
</td>
<td style="text-align:right;">
0.6780636
</td>
<td style="text-align:right;">
0.8234462
</td>
<td style="text-align:right;">
-5.8854
</td>
<td style="text-align:right;">
-2.6575
</td>
<td style="text-align:right;">
0.0000025
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
167
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
X4-X5
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.1345843
</td>
<td style="text-align:right;">
0.7007521
</td>
<td style="text-align:right;">
0.8371094
</td>
<td style="text-align:right;">
-1.5061
</td>
<td style="text-align:right;">
1.7753
</td>
<td style="text-align:right;">
0.8830556
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
167
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Interaction
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.4549679
</td>
<td style="text-align:right;">
0.5524247
</td>
<td style="text-align:right;">
0.7432528
</td>
<td style="text-align:right;">
-1.9117
</td>
<td style="text-align:right;">
1.0018
</td>
<td style="text-align:right;">
0.5976862
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
166
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
X4
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-3.3168676
</td>
<td style="text-align:right;">
0.4551012
</td>
<td style="text-align:right;">
0.6746119
</td>
<td style="text-align:right;">
-4.6391
</td>
<td style="text-align:right;">
-1.9947
</td>
<td style="text-align:right;">
0.0000538
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
166
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
X5
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-3.9151131
</td>
<td style="text-align:right;">
0.4671708
</td>
<td style="text-align:right;">
0.6834989
</td>
<td style="text-align:right;">
-5.2547
</td>
<td style="text-align:right;">
-2.5755
</td>
<td style="text-align:right;">
0.0000022
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
166
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
X4-X5
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-0.1432777
</td>
<td style="text-align:right;">
0.5245870
</td>
<td style="text-align:right;">
0.7242838
</td>
<td style="text-align:right;">
-1.5628
</td>
<td style="text-align:right;">
1.2763
</td>
<td style="text-align:right;">
0.8663046
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
166
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Interaction
</td>
</tr>
</tbody>
</table>

Here we see that the interaction between X4 and X5 was consistently
found to have the highest synergistic interaction across the folds.
Therefore, for our pooled parameter var 1 represents the pooled effects
of shifting X4, var 2 represents the pooled effects of shifting X5,
joint is X4 and X5 together and the interaction represents the
interaction effect for these two variables.

Next we’ll look at the k-fold antagonistic interactions:

``` r
k_fold_antagonism_effects$`Rank 1` %>%
  kbl(caption = "K-fold Antagonistic Stochastic Intervention Results") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
K-fold Antagonistic Stochastic Intervention Results
</caption>
<thead>
<tr>
<th style="text-align:right;">
Rank
</th>
<th style="text-align:right;">
Psi
</th>
<th style="text-align:right;">
Variance
</th>
<th style="text-align:right;">
SE
</th>
<th style="text-align:right;">
Lower CI
</th>
<th style="text-align:right;">
Upper CI
</th>
<th style="text-align:right;">
P-value
</th>
<th style="text-align:right;">
Fold
</th>
<th style="text-align:right;">
N
</th>
<th style="text-align:right;">
Delta Exposure 1
</th>
<th style="text-align:right;">
Delta Exposure 2
</th>
<th style="text-align:left;">
Type
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
13.542234
</td>
<td style="text-align:right;">
0.4842421
</td>
<td style="text-align:right;">
0.6958751
</td>
<td style="text-align:right;">
12.1783
</td>
<td style="text-align:right;">
14.9061
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
167
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
X1
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
3.487427
</td>
<td style="text-align:right;">
0.7311449
</td>
<td style="text-align:right;">
0.8550701
</td>
<td style="text-align:right;">
1.8115
</td>
<td style="text-align:right;">
5.1633
</td>
<td style="text-align:right;">
0.0001623
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
167
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
X7
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
10.366036
</td>
<td style="text-align:right;">
0.3799426
</td>
<td style="text-align:right;">
0.6163949
</td>
<td style="text-align:right;">
9.1579
</td>
<td style="text-align:right;">
11.5741
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
167
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
X1-X7
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-6.663625
</td>
<td style="text-align:right;">
0.6288383
</td>
<td style="text-align:right;">
0.7929932
</td>
<td style="text-align:right;">
-8.2179
</td>
<td style="text-align:right;">
-5.1094
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
167
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Interaction
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-1.245423
</td>
<td style="text-align:right;">
1.5739266
</td>
<td style="text-align:right;">
1.2545623
</td>
<td style="text-align:right;">
-3.7043
</td>
<td style="text-align:right;">
1.2135
</td>
<td style="text-align:right;">
0.2661755
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
167
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
X1
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
2.963629
</td>
<td style="text-align:right;">
0.8560357
</td>
<td style="text-align:right;">
0.9252220
</td>
<td style="text-align:right;">
1.1502
</td>
<td style="text-align:right;">
4.7770
</td>
<td style="text-align:right;">
0.0020626
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
167
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
X7
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
7.220341
</td>
<td style="text-align:right;">
0.5888704
</td>
<td style="text-align:right;">
0.7673789
</td>
<td style="text-align:right;">
5.7163
</td>
<td style="text-align:right;">
8.7244
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
167
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
X1-X7
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
5.502135
</td>
<td style="text-align:right;">
2.0114489
</td>
<td style="text-align:right;">
1.4182556
</td>
<td style="text-align:right;">
2.7224
</td>
<td style="text-align:right;">
8.2819
</td>
<td style="text-align:right;">
0.0000038
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
167
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Interaction
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
16.279113
</td>
<td style="text-align:right;">
0.4910244
</td>
<td style="text-align:right;">
0.7007314
</td>
<td style="text-align:right;">
14.9057
</td>
<td style="text-align:right;">
17.6525
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
166
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
X1
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
2.878816
</td>
<td style="text-align:right;">
0.4632692
</td>
<td style="text-align:right;">
0.6806388
</td>
<td style="text-align:right;">
1.5448
</td>
<td style="text-align:right;">
4.2128
</td>
<td style="text-align:right;">
0.0004840
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
166
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
X7
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
7.522615
</td>
<td style="text-align:right;">
0.3101366
</td>
<td style="text-align:right;">
0.5568991
</td>
<td style="text-align:right;">
6.4311
</td>
<td style="text-align:right;">
8.6141
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
166
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
X1-X7
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-11.635314
</td>
<td style="text-align:right;">
0.4261335
</td>
<td style="text-align:right;">
0.6527890
</td>
<td style="text-align:right;">
-12.9148
</td>
<td style="text-align:right;">
-10.3559
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
166
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Interaction
</td>
</tr>
</tbody>
</table>

Here, we see that in all the folds the X1-X7 has the strongest
antagonistic relationship. X1-X7 was found in all the folds and
therefore the oracle parameter is interpreted the same as we found in
the synergy results. Which is here:

``` r
pooled_antagonism_effects$`Rank 1` %>%
  kbl(caption = "Rank 1 Antagonism Stochastic Intervention Results") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Rank 1 Antagonism Stochastic Intervention Results
</caption>
<thead>
<tr>
<th style="text-align:left;">
Rank
</th>
<th style="text-align:right;">
Psi
</th>
<th style="text-align:right;">
Variance
</th>
<th style="text-align:right;">
SE
</th>
<th style="text-align:right;">
Lower CI
</th>
<th style="text-align:right;">
Upper CI
</th>
<th style="text-align:right;">
P-value
</th>
<th style="text-align:left;">
Fold
</th>
<th style="text-align:right;">
N
</th>
<th style="text-align:right;">
Delta Exposure 1
</th>
<th style="text-align:right;">
Delta Exposure 2
</th>
<th style="text-align:left;">
Type
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Rank 1
</td>
<td style="text-align:right;">
13.215850
</td>
<td style="text-align:right;">
0.2422382
</td>
<td style="text-align:right;">
0.4921770
</td>
<td style="text-align:right;">
12.2512
</td>
<td style="text-align:right;">
14.1805
</td>
<td style="text-align:right;">
0.0e+00
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:right;">
500
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Var 1
</td>
</tr>
<tr>
<td style="text-align:left;">
Rank 1
</td>
<td style="text-align:right;">
3.253814
</td>
<td style="text-align:right;">
0.2683451
</td>
<td style="text-align:right;">
0.5180204
</td>
<td style="text-align:right;">
2.2385
</td>
<td style="text-align:right;">
4.2691
</td>
<td style="text-align:right;">
6.2e-06
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:right;">
500
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Var 2
</td>
</tr>
<tr>
<td style="text-align:left;">
Rank 1
</td>
<td style="text-align:right;">
9.138029
</td>
<td style="text-align:right;">
0.1578888
</td>
<td style="text-align:right;">
0.3973522
</td>
<td style="text-align:right;">
8.3592
</td>
<td style="text-align:right;">
9.9168
</td>
<td style="text-align:right;">
0.0e+00
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:right;">
500
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Joint
</td>
</tr>
<tr>
<td style="text-align:left;">
Rank 1
</td>
<td style="text-align:right;">
-7.331635
</td>
<td style="text-align:right;">
0.3119851
</td>
<td style="text-align:right;">
0.5585563
</td>
<td style="text-align:right;">
-8.4264
</td>
<td style="text-align:right;">
-6.2369
</td>
<td style="text-align:right;">
0.0e+00
</td>
<td style="text-align:left;">
Pooled TMLE
</td>
<td style="text-align:right;">
500
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Interaction
</td>
</tr>
</tbody>
</table>

So we see the interaction effect is negative, -7.33, and represents the
pooled interaction effects across the folds which are all X1-X7.

Overall, this package provides implementation of estimation a
non-parametric definition of interaction. We define positive values as
synergy meaning the expected outcome under joint shift is much larger
compared to individual addivitive effects. Likewise, we define
antagonism as negative effects, the joint value being lower than the
additive effects.

In this NIEHS data set w correctly identify the strongest individual
effects in positive and negative directions and identify exposure
relationships consistently for our definition of synergy and antagonism.

------------------------------------------------------------------------

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/blind-contours/InterXshift/issues).
Further details on filing issues are provided in our [contribution
guidelines](https://github.com/blind-contours/%20InterXshift/main/contributing.md).

------------------------------------------------------------------------

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/blind-contours/InterXshift/blob/master/CONTRIBUTING.md)
prior to submitting a pull request.

------------------------------------------------------------------------

## Citation

After using the `InterXshift` R package, please cite the following:

------------------------------------------------------------------------

## Related

- [R/`tmle3shift`](https://github.com/tlverse/tmle3shift) - An R package
  providing an independent implementation of the same core routines for
  the TML estimation procedure and statistical methodology as is made
  available here, through reliance on a unified interface for Targeted
  Learning provided by the [`tmle3`](https://github.com/tlverse/tmle3)
  engine of the [`tlverse` ecosystem](https://github.com/tlverse).

- [R/`medshift`](https://github.com/nhejazi/medshift) - An R package
  providing facilities to estimate the causal effect of stochastic
  treatment regimes in the mediation setting, including classical (IPW)
  and augmented double robust (one-step) estimators. This is an
  implementation of the methodology explored by Dı́az and Hejazi (2020).

- [R/`haldensify`](https://github.com/nhejazi/haldensify) - A minimal
  package for estimating the conditional density treatment mechanism
  component of this parameter based on using the [highly adaptive
  lasso](https://github.com/tlverse/hal9001) (Coyle, Hejazi, Phillips,
  et al. 2022; Hejazi, Coyle, and van der Laan 2020) in combination with
  a pooled hazard regression. This package implements a variant of the
  approach advocated by Dı́az and van der Laan (2011).

------------------------------------------------------------------------

## Funding

The development of this software was supported in part through grants
from the

------------------------------------------------------------------------

## License

© 2020-2022 [David B. McCoy](https://davidmccoy.org)

The contents of this repository are distributed under the MIT license.
See below for details:

    MIT License
    Copyright (c) 2020-2022 David B. McCoy
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

------------------------------------------------------------------------

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-coyle-sl3-rpkg" class="csl-entry">

Coyle, Jeremy R, Nima S Hejazi, Ivana Malenica, Rachael V Phillips, and
Oleg Sofrygin. 2022. “<span class="nocase">sl3</span>: Modern Machine
Learning Pipelines for Super Learning.”
<https://doi.org/10.5281/zenodo.1342293>.

</div>

<div id="ref-coyle-hal9001-rpkg" class="csl-entry">

Coyle, Jeremy R, Nima S Hejazi, Rachael V Phillips, Lars W van der Laan,
and Mark J van der Laan. 2022. “<span class="nocase">hal9001</span>: The
Scalable Highly Adaptive Lasso.”
<https://doi.org/10.5281/zenodo.3558313>.

</div>

<div id="ref-diaz2020causal" class="csl-entry">

Dı́az, Iván, and Nima S Hejazi. 2020. “Causal Mediation Analysis for
Stochastic Interventions.” *Journal of the Royal Statistical Society:
Series B (Statistical Methodology)* 82 (3): 661–83.
<https://doi.org/10.1111/rssb.12362>.

</div>

<div id="ref-diaz2011super" class="csl-entry">

Dı́az, Iván, and Mark J van der Laan. 2011. “Super Learner Based
Conditional Density Estimation with Application to Marginal Structural
Models.” *The International Journal of Biostatistics* 7 (1): 1–20.

</div>

<div id="ref-diaz2012population" class="csl-entry">

———. 2012. “Population Intervention Causal Effects Based on Stochastic
Interventions.” *Biometrics* 68 (2): 541–49.

</div>

<div id="ref-haneuse2013estimation" class="csl-entry">

Haneuse, Sebastian, and Andrea Rotnitzky. 2013. “Estimation of the
Effect of Interventions That Modify the Received Treatment.” *Statistics
in Medicine* 32 (30): 5260–77.

</div>

<div id="ref-hejazi2020hal9001-joss" class="csl-entry">

Hejazi, Nima S, Jeremy R Coyle, and Mark J van der Laan. 2020.
“<span class="nocase">hal9001</span>: Scalable Highly Adaptive Lasso
Regression in R.” *Journal of Open Source Software* 5 (53): 2526.
<https://doi.org/10.21105/joss.02526>.

</div>

<div id="ref-mccoy2023semiparametric" class="csl-entry">

McCoy, David B., Alan E. Hubbard, Alejandro Schuler, and Mark J. van der
Laan. 2023. “Semi-Parametric Identification and Estimation of
Interaction and Effect Modification in Mixed Exposures Using Stochastic
Interventions.” <https://arxiv.org/abs/2305.01849>.

</div>

</div>
