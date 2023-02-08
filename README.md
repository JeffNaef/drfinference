# Condifence and Uncertainty Assessment for Distributional Random Forests: Code

This repo contains the code for generating the figures in the paper "Condifence and Uncertainty Assessment for Distributional Random Forests" by Jeffrey N&auml;f, Corinne Emmenegger, Peter B&uuml;hlmann, and Nicolai Meinshausen. 

## Target parameter estimation

Run the file `target-param.R` and specify the following variables:

- Figure 2 (CATE with homogeneous treatment effect and observed confounding): `do_copula = FALSE`, `do_quantile = FALSE`, `sim = 1`
- Figure 3 (CATE with heterogeneous treatment effect and observed confounding): `do_copula = FALSE`, `do_quantile = FALSE`, `sim = 3`
- Figure 4 and 5 (conditional quantiles): `do_copula = FALSE`, `do_quantile = TRUE`
- Figure 6 (conditional correlation): `do_copula = TRUE`, `do_quantile = FALSE`

## Witness function for conditional distributional treatment effect

- Figure 7 (investigate where the treatment and control distributions differ with data generating mechanism from Figure 2): run the file `distr-difference.R` with `sim = 1`. To get the accompanying information for this figure for the data generating mechanism as in Figure 3, set `sim = 3`.
- Figure 8 (simultaneous confidence bands for conditional witness function): run the file `witnessfunc.R` with `sim = 1` to obtain Figure 8 (a) and with `sim = 3` to obtain Figure 8 (b). 

The results after running any of the files will be saved in the folder `Results`.