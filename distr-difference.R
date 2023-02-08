#' based on run_witnessfunc_v2_sim1_runner.R and run_witnessfunc_v2_sim3_runner.R
#' data-foo-codite: based on dataGen_codite.R
#' plot-codite: based on plot_codite_v2.R
#' drf-foo: based on drfnew_v2.R

# load libraries and files
library(kernlab)
library(drf)
library(Matrix)
library(dplyr)
library(doParallel)
library(doRNG)
library(parallel)
library(foreach)
library(ggplot2)
library(ggpubr)
print(sessionInfo())

source("drf-foo.R")
source("data-foo-codite.R")
source("plot-codite.R")

### random seed and save current working directory
RNGkind("L'Ecuyer-CMRG")
set.seed(1)
old_wd <- getwd()

### initialize parameters
num.trees <- 1000 # number of trees of each individual halfsampling forest
B <- 300 # number of halfsampling DRF's (based on these, variance is computed)
reps <- 20  # number of repetitions for simulation
n_vec <- 10000
save_filename <- "res"
sampling <- "binomial" # sampling mechanism to generate halfsamples
predictvar <- TRUE
cores <- 20
no_splitup <- 50 # total repetitions = no_splitup * reps and for
                 # each set of reps many repetitions, use cores many cores
                 # in parallelization (no_splitup >= 2 required)

# specification of DRF
num.features <- 10
honesty <- TRUE
compute.oob.predictions <- FALSE
response.scaling <- TRUE
min.node.size <- 5
complement_half_sampling <- TRUE
beta <- 0.95 # fraction of observations used to fit one tree
alpha_ks <- 0.05 # level for ks test of nulldist0 = nulldist1

# cate
sim <- 3 # data generation. either 1 (no treatment effect) or 3 (treatment effect)
p <- 5
q <- 1
d <- 2
data_foo <- data_foo_cate
get_truth <- get_truth_cate

# save setting
save_setting(num.trees = num.trees, beta = beta,
             B = B, predictvar = predictvar, sim = sim,
             save_filename = save_filename,
             subdir = "Results",
             reps = reps,
             num.features = num.features,
             honesty = honesty,
             response.scaling = response.scaling,
             compute.oob.predictions = compute.oob.predictions,
             min.node.size = min.node.size,
             sampling = sampling,
             data_foo = data_foo,
             n_vec = n_vec,
             alpha_ks = alpha_ks,
             d = d, q = q, p = p,
             get_truth = get_truth)

for (kk in seq_len(no_splitup)) {
  print(paste0("kk = ", kk))

  # initiate parallelization
  registerDoParallel(cores = cores)
  seeds <- 1:reps + (kk - 1) * reps
  (time_start <- Sys.time()) # save start running time

  # do parallelized simultions
  for (n in n_vec) {
    print(paste0("n = ", n, ", ", Sys.time()))
    s <- n ^ beta
    res <- foreach(rep = 1:reps, .combine = rbind, .packages = c("drf", "Matrix", "kernlab", "dplyr")) %dopar% {

      set.seed(seeds[rep])

      if ((rep / reps * 100) %% 10 == 0) {
        cat(paste0(rep / reps * 100, "%  "))
      }

      # simulate data
      data <- data_foo(sim = sim, n = n, p = p, q = q, d = d)

      # Fit two DRFs for the two classes
      DRF0 <-
        drfCI(X = data$X[data$W == 0, , drop = FALSE],
              Y = data$Y[data$W == 0, , drop = FALSE],
              B = B, num.trees = num.trees)
      DRF1 <-
        drfCI(X = data$X[data$W == 1, , drop = FALSE],
              Y = data$Y[data$W == 1, , drop = FALSE],
              B = B, num.trees = num.trees)
      # predict the DRFs on testdata x
      DRFpred0 <- predictdrf(DRF0, x = data$x)
      DRFpred1 <- predictdrf(DRF1, x = data$x)

      # kernel
      bandwidth_Y <- drf:::medianHeuristic(data$Y)
      k_Y <- rbfdot(sigma = bandwidth_Y)

      K1 <- kernelMatrix(k_Y, data$Y[data$W == 1], y = data$Y[data$W == 1])
      K0 <- kernelMatrix(k_Y, data$Y[data$W == 0], y = data$Y[data$W == 0])
      K <- kernelMatrix(k_Y, data$Y[data$W == 0], y = data$Y[data$W == 1])

      # ||mu_{n,1} - mu_{n,0}||^2
      teststat <-
        diag(DRFpred0$weights %*% tcrossprod(K0, DRFpred0$weights)) +
        diag(DRFpred1$weights %*% tcrossprod(K1, DRFpred1$weights)) -
        2 * diag(DRFpred0$weights %*% tcrossprod(K,  DRFpred1$weights)) %>%
        as.numeric()

      # simulated null distribution of test statistic
      nulldist <- sapply(seq_len(length(DRFpred1$weightsb)), function(j) {
        # iterate over class 1

        wb0 <- DRFpred0$weightsb[[j]]
        wb1 <- DRFpred1$weightsb[[j]]

        diag((wb0 - DRFpred0$weights) %*% tcrossprod(K0, wb0 - DRFpred0$weights)) +
          diag((wb1 - DRFpred1$weights) %*% tcrossprod(K1, wb1 - DRFpred1$weights)) -
          2 * diag((wb0 - DRFpred0$weights) %*% tcrossprod(K, (wb1 - DRFpred1$weights))) %>%
          as.numeric()

      })
      one_dim_x <- dim(data$x)[1] == 1
      if (one_dim_x) {
        nulldist <- matrix(nulldist, ncol = B)
      }

      tibble(pvalue = sapply(seq_len(length(teststat)),
                             function(i) {
                               (sum(nulldist[i, ] > teststat[i]) + 1) / (length(nulldist[i, ]) + 1)
                             }),
             teststat = teststat,
             x = data$x[, 1] %>%
               as.numeric()) %>%
        return()

    }

    print(res)
    save(res, file = paste0(save_filename, "_n_", n, "_run_", kk, ".RData"))
  }
}

do_plots(runout = no_splitup)