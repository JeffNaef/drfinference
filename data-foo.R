# when running a simulation, this function saves all parameter choices
save_setting <- function(num.trees, beta, l = NULL,
                         B, predictvar, sim, dataset, save_filename,
                         subdir = "Results",
                         reps, x,
                         num.features,
                         honesty,
                         response.scaling,
                         compute.oob.predictions,
                         min.node.size,
                         param_foo,
                         thetazero,
                         sampling,
                         data_foo,
                         get_truth,
                         n_vec,
                         do_CF, do_quantile,
                         d, p, q, probs,
                         SC, do_copula, alpha) {
  now <- Sys.time()
  print(now)
  now.date <- format(now, "%d_%b_%Y")
  now.time <- format(now, "%H_%M_%S")
  dirname <- paste(subdir, "/pics__", now.date, "__", now.time, sep = "")

  if (!dir.exists(subdir)) {
    dir.create(subdir)
  }
  if (!dir.exists(dirname)) {
    dir.create(dirname)
  }
  setwd(dirname)

  setup <- list(d = d, p = p, q = q,
                num.trees = num.trees, beta = beta,
                B = B, l = l, predictvar = predictvar, sim = sim,
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
                do_CF = do_CF, do_quantile = do_quantile,
                probs = probs,
                do_center = do_center,
                do_copula = do_copula,
                SC = SC,
                alpha = alpha)
  save(setup, file = "preliminary.RData")
  dput(setup, file = "preliminary.txt")
}

# data simulation scenarios from DRF paper from Section D.2.1 for cate
fooY <- function(x) {
  1 + (1 + exp(-20 * (x - 1 / 3))) ^ (-1)
}
expit <- function(x) {
  exp(x) / (1 + exp(x))
}
data_foo_cate <- function(sim, n, q, p, d) {
  #' p = dimension of X
  #' d = dimension of Y
  #' q = dimension of theta

  X <- matrix(runif(p * n), ncol = p, nrow = n, byrow = T)

  if (sim == 1) { # First simulation (confounder affecting W and Y)
    stopifnot(q == 1)
    stopifnot(d == 2)
    stopifnot(p >= 3)
    W <- rbinom(n, 1, (1 + dbeta(X[, 3], 2, 4)) / 4)
    Y <-
      cbind(2 * (X[, 3] - 0.5) + rnorm(n),
            W)

  } else if (sim == 2) { # Second simulation (heterogenous treatment effect)
    stopifnot(d == 2)
    stopifnot((q == 1) | (q == 2))
    W <- rbinom(n, 1, 0.5)
    etaX <- apply(X, 2, function(x) fooY(x))
    Y <-
      cbind((W - 1 / 2) * etaX[, 1] * etaX[, 2] + rnorm(n),
            W)

  } else if (sim == 3) { # Third simulation (heterogeneous treatment effect + confounding)
    stopifnot(d == 2)
    stopifnot((q == 1) | (q == 2))
    stopifnot(p >= 3)
    Y <- matrix(rnorm(d * n), ncol = d, nrow = n, byrow = T)
    px3 <- 1 / 4 * (1 + dbeta(X[, 3], shape1 = 2, shape2 = 4))
    etaX <- apply(X, 2, function(x) fooY(x))
    # Y[, 2] is the treatment assignment W:
    Y[, 2] <- rbinom(n = n, size = 1, prob = px3)
    Y[, 1] <- 2 * (X[, 3] - 1 / 2) - 1 / 2 * etaX[, 1] * etaX[, 2] +
      Y[, 2] * etaX[, 1] * etaX[, 2] + rnorm(n)

  } else if (sim == 4) {
    stopifnot(d == 2)
    stopifnot(p >= 2)
    stopifnot((q == 1) | (q == 2))
    W <- rbinom(n = n, size = 1, prob = expit(4 * X[, 2] - 2))
    Y <-
      cbind(100 * X[, 2] ^ 2 + (W - 0.5) * sin(3 * X[, 1]) + rnorm(n),
            W)
  }

  list(X = X, Y = Y, W = NULL)
}
# function to get thetazero and the conditional x with cate
get_truth_cate <- function(sim, q, p) {

  if (sim == 1) {
    x <- matrix(c(c(0.7, 0.3, 0.5), runif(p))[seq_len(p)], ncol = p)
    thetazero <- 2 * (x[3] - 0.5)

  } else if (sim == 2) {
    x <- matrix(c(c(0.7, 0.3, 0.5), runif(p))[seq_len(p)], ncol = p)
    etax <- apply(x, 2, function(x) fooY(x))
    thetazero <- etax[1] * etax[2]

    if (q == 2) { # for 2-dimensional theta, consider intercept and slope
      thetazero <-
        c(-1 / 2 * thetazero, thetazero)
    }

  } else if (sim == 3) {
    x <- matrix(c(c(0.7, 0.3, 0.5), runif(p))[seq_len(p)], ncol = p)
    etax <- apply(x, 2, function(x) fooY(x))
    thetazero <- etax[1] * etax[2]

    if (q == 2) {
      thetazero <- c(2 * (x[, 3] - 1 / 2) - 1 / 2 * etax[1] * etax[2], thetazero)
    }

  } else if (sim == 4) {
    x <- matrix(c(c(0.7, 0.3, 0.5), runif(p))[seq_len(p)], ncol = p)
    thetazero <- sin(3 * x[1])

    if (q == 2) {
      thetazero <- c(100 * x[2] ^ 2 - 0.5 * sin(3 * x[1]))
    }
  }

  list(x = x, thetazero = thetazero)
}

# data for quantile example (as in the DRF paper)
data_foo_quant <- function(dataset = "synthetic2", n = 5000, p = 39,
                           meanShift = 0.8, sdShift = 1, probs = c(0.3, 0.5, 0.7)) {
  #' data generating mechanisms corresponding to quantile analysis in DRF paper
  #' corresponding to figure 4
  if (dataset == "synthetic1") {
    x <- runif(n, -1, 1)
    y <- rnorm(n, mean = meanShift * (x > 0))
    X <- cbind(x, matrix(runif(n * p), ncol = p))

  } else if (dataset == "synthetic2") {
    x <- runif(n, -1, 1)
    y <- rnorm(n, sd = 1 + sdShift * (x > 0))
    X <- cbind(x, matrix(runif(n * p), ncol = p))

  } else if (dataset == "synthetic3") {
    x <- runif(n, -1, 1)
    y <- ifelse(x >= 0, rexp(n = n, 1), rnorm(n, 1, 1))
    X <- cbind(x, matrix(runif(n * p), ncol = p))

  } else {
    stop("no valid dataset name in function genData")
  }

  list(Y = as.matrix(y, ncol = 1), X = X)
}
# get true thetazero and conditional x
get_truth_quant <- function(dataset, p = 39, probs = c(0.3, 0.5, 0.7),
                            meanShift = 0.8, sdShift = 1) {
  xtest <- -10:10 / 51
  Xtest <- matrix(0, nrow = length(xtest), ncol = p)
  thetazero <- if (dataset == "synthetic1") {
    sapply(seq_len(length(xtest)),
           function(i) qnorm(probs, mean = meanShift * (xtest[i] > 0))) %>%
      t()
  } else if (dataset == "synthetic2") {
    sapply(seq_len(length(xtest)),
           function(i) qnorm(probs, sd = 1 + sdShift * (xtest[i] > 0))) %>%
      t()
  } else if (dataset == "synthetic3") {
    sapply(seq_len(length(xtest)),
           function(i) if (xtest[i] >= 0) {
             qexp(probs, 1)
           } else {
             qnorm(probs, 1, 1)
           }) %>%
      t()
  }

  list(x = cbind(xtest, Xtest), thetazero = thetazero)
}

# data for copula example: estimating conditional correlation between Y1 and Y2
# according to DRF paper
# SC = 0: normal copula example with N(0,1) marginals, change of correlation of Y1 and Y2 depends on X1
# SC = 1: t-copula, X1 is tail parameter, X2 is covariance and X3 is marginal
# SC = 2: normal copula example with N(0,1) marginals, toeplitz covariance rho depends on X1
# SC = 3: normal copula example with N(0,1) marginals, equicorrelation covariance rho depends on X1
gen_copula <- function(xx, SC, d) { # helper function
  if (SC == 0) {
    # copula
    normCop <- normalCopula(param = c(xx[1]), dim = 2)
    # margins
    paramMargins <- list(list(mean = 0, sd = 1), list(mean = 0, sd = 1))
    mdvNorm <- mvdc(copula = normCop, margins = c("norm", "norm"), paramMargins = paramMargins)
    # gen
    c(rMvdc(n = 1, mvdc = mdvNorm), rnorm(d - 2))
  }
  else if (SC == 1) {
    # copula
    tCop <- tCopula(xx[2],  df = ifelse(xx[1] <= (-1 + 2 / 3), 1, ifelse(xx[1] <= (-1 + 4 / 3), 3, 10)))
    # margins
    if (xx[3] <= 0) {
      margins <- c("norm", "norm")
      paramMargins <- list(list(mean = 0, sd = 1), list(mean = 0, sd = 1))
    } else {
      margins <- c("norm", "exp")
      paramMargins <- list(list(mean = 0, sd = 1), list(rate = 1))
    }
    mdvT <- mvdc(copula = tCop, margins = margins, paramMargins = paramMargins)
    # gen
    c(rMvdc(n = 1, mvdc = mdvT), rnorm(d - 2))

  }
  else if (SC == 2) {
    mu <- rep(0, d)
    rho <- xx[1]
    Sigma <- toeplitz(rho ^ (pmin(0:(d - 1), d - 0:(d - 1))))
    mvrnorm(mu = mu, Sigma = Sigma)
  }
  else if (SC == 3) {
    mu <- rep(0, d)
    rho <- xx[1]
    Sigma <- toeplitz(c(1, rep(rho, d - 1)))
    mvrnorm(mu = mu, Sigma = Sigma)
  }
}
# actual data generating fnction
data_foo_copula <- function(SC = 3, n = 5000, q = 1, p = 30, d = 5) {
  stopifnot(q == 1)
  # GENERATE DATA

  # predictors
  X <- if (SC == 3 || SC == 2) {
    matrix(runif(n * p, min = -1, max = 1), ncol = p)
  } else if (SC == 0 || SC == 1) {
    matrix(runif(n * p, min = -1, max = 1), ncol = p)
  }
  # responses
  Y <- t(apply(X, 1, gen_copula, SC = SC, d = d))
  colnames(Y) <- paste0("Y", seq_len(d))

  # return
  list(Y = Y, X = X)
}
# get true thetazero and conditioning x for copula
get_truth_copula <- function(SC, p = 30) {
  a <- seq(0.000108375, 0.999698272, length.out = 10)
  x <- matrix(sort(c(-a, a)), ncol = 1)
  x_test <- cbind(x, matrix(0.5005142, nrow = nrow(x), ncol = p - 1))
  thetazero <- if (SC == 0 | SC == 1) {
    stop("not implemented")
  } else if (SC == 2 | SC == 3) {
    x
  }

  list(x = x_test, thetazero = thetazero)
}
