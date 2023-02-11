#' this file corresponds to runboot.R
#' data-foo corresponds to dataGen
#' helper-foo corresponds to helpers
#' plot-target-param corresponds to plot_v2S1_v3

### load files and packages
source("helper-foo.R")
source("data-foo.R")
source("plot-target-param.R")
library(doParallel)
library(doRNG)
library(parallel)
library(Matrix)
library(kernlab)
library(drf)
library(grf)
library(tidyverse)
library(expm)
library(mixtools)
library(foreach)
library(doParallel)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(causalToolbox)
library(locfit)
library(ranger)
library(copula)
library(MASS)
print(sessionInfo())

### random seed and save current working directory
RNGkind("L'Ecuyer-CMRG")
set.seed(1)
old_wd <- getwd()

### choose what to estimate 
### if copula and quantile are both false, cate is done
do_copula <- FALSE
do_quantile <- FALSE

# cate
sim <- 3 # data generation. To compute cate, choose 1 (homogeneous treatment effect and observed confounding)
         # or choose 3 (heterogeneous treatment effect and observed confounding)
do_CF <- TRUE # compute GRF 
if (do_quantile | do_copula) {
  do_CF <- FALSE
}

# copula: conditional correlation cor(Y1, Y2 | X)
SC <- 3 # data generation

# quantile: conditional quantile
dataset <- "synthetic1" # data generation
probs <- c(0.1, 0.5, 0.9) # which conditional quantiles


### initialize parameters
p <- if (do_quantile) { # dimension of X
  4
} else if (do_copula) {
  5
} else {
  5
}
d <- if (do_quantile) { # dimension of Y
  1
} else if (do_copula) {
  2
} else {
  2
}
q <- if (do_quantile) { # dimension of theta
  length(probs) #3
} else if (do_copula) {
  1
} else {
  1
}
num.trees <- 1000 # number of trees of each individual halfsampling forest
B <- 100 # number of halfsampling DRF's (based on these, variance is computed)
reps <- 1000 # number of repetitions for simulation
n_vec <- if (do_quantile) { # which sample sizes are considered
  5000
} else if (do_copula) {
  c(312, 625, 1250, 2500, 5000)
} else {
  c(312, 625, 1250, 2500, 5000)
}
alpha <- 0.05 # level for testing
sampling <- "binomial" # sampling mechanism to generate halfsamples
predictvar <- TRUE # predict the variance

save_filename <- "res" # prename used to save files
cores <- 84 # number of cores used for parallel computing

# specification of DRF
num.features <- 10
honesty <- TRUE
compute.oob.predictions <- FALSE
response.scaling <- FALSE
min.node.size <- 5
complement_half_sampling <- TRUE
beta <- 0.9 # fraction of observations used to fit one tree


### get data generating function
data_foo <- if (do_quantile) {
  data_foo_quant
} else if (do_copula) {
  data_foo_copula
} else {
  data_foo_cate
}
### get true values function
get_truth <- if (do_quantile) {
  get_truth_quant
} else if (do_copula) {
  get_truth_copula
} else {
  get_truth_cate
}
# get parameter computation function
param_foo <- if (do_quantile) { # conditional quantile
  function(weights, X, Y, W = NULL, x, probs = c(0.1, 0.5, 0.9)) {
    drf:::weighted.quantile(c(Y), weights, probs = probs)
  }
} else if (do_copula) { # conditional correlation
  function(weights, X, Y, W = NULL, x, probs = NULL) {
    # cor(Y1, Y2 | X = x)
    cov.wt(Y[, 1:2], wt = weights %>% as.vector(), cor = TRUE)$cor[1, 2] # from library stats
  }
} else if (q == 1) { # cate (slope only)
  function(weights, X, Y, W = NULL, x, probs = NULL) {
    unname(coef(lm(Y[, 1] ~ Y[, 2:ncol(Y)], weights = as.numeric(weights)))[-1])
  }
} else if (q == 2) { # cate (slope and intercept)
  function(weights, X, Y, W = NULL, x, probs = NULL) {
    unname(coef(lm(Y[, 1] ~ Y[, 2:ncol(Y)], weights = as.numeric(weights))))
  }
}

# extract thetazero and x
tmp <- if (do_quantile) {
  get_truth(dataset = dataset, p = p, probs = probs)
} else if (do_copula) {
  get_truth(SC = SC, p = p)
} else {
  get_truth(sim = sim, q = q, p = p)
}
x <- tmp$x
thetazero <- tmp$thetazero

# save the current setting
save_setting(d = d, p = p, q = q,
             num.trees = num.trees, beta = beta,
             B = B,
             predictvar = predictvar, sim = sim,
             dataset = dataset,
             save_filename = save_filename,
             reps = reps, x = x,
             num.features = num.features,
             honesty = honesty,
             response.scaling = response.scaling,
             compute.oob.predictions = compute.oob.predictions,
             min.node.size = min.node.size,
             param_foo = param_foo,
             thetazero = thetazero,
             sampling = sampling,
             data_foo = data_foo,
             get_truth = get_truth,
             n_vec = n_vec,
             do_CF = do_CF, 
             do_quantile = do_quantile,
             probs = probs,
             do_copula = do_copula,
             SC = SC,
             alpha = alpha)

# initiate parallelization
registerDoParallel(cores = cores)
seeds <- 1:reps
(time_start <- Sys.time()) # save start running time

# do parallelized simultions
for (n in n_vec) {
  print(paste0("n = ", n, ", ", Sys.time()))
  s <- n ^ beta
  res <- foreach(rep = 1:reps, .combine = rbind, .packages = c("drf", "Matrix")) %dopar% { #.export=ls(.GlobalEnv),
    
    set.seed(seeds[rep])
    
    if ((rep / reps * 100) %% 10 == 0) {
      cat(paste0(rep / reps * 100, "%  "))
    }
    
    # simulate data
    data <- if (do_quantile) {
      data_foo(n = n, probs = probs, dataset = dataset, p = p)
    } else if (do_copula) {
      data_foo(SC = SC, n = n, q = q, p = p, d = d)
    } else {
      data_foo(sim, n = n, d = d, q = q, p = p)
    }
    
    # DRF
    tmp <-
      bootvar(Y = data$Y, X = data$X, param_foo, x = x,
              num.features = num.features,
              honesty = honesty,
              num.trees = num.trees, B = B,
              beta = beta,
              response.scaling = response.scaling,
              compute.oob.predictions = compute.oob.predictions,
              min.node.size = min.node.size,
              sampling = sampling)
    # extract estimator and variance
    to_return <- list(thetahat = tmp$thetahat, var_est = tmp$var_est)
    
    # causal forest
    if (do_CF) {
      CF_fit <- causal_forest(X = data$X, Y = data$Y[, 1], W = data$Y[, 2], 
                              num.trees = min(B * num.trees, 50000), 
                              sample.fraction = 1 - beta)
      CF_predict <- predict(CF_fit, x, estimate.variance = TRUE) #computes CATE
      CATE_CF <- CF_predict[, "predictions"]
      CF_var <- CF_predict[, "variance.estimates"]
      
      to_return <-
        c(to_return, list(thetahat_CF = CATE_CF, var_est_CF = CF_var))
    }
    
    return(to_return)
  }
  
  save(res, file = paste0(save_filename, "_n_", n, ".RData"))
}

# print computing time
print(Sys.time() - time_start)
# make plots
do_plots(do_copula = do_copula)
# set initial wd
setwd(old_wd)
