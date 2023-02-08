# when running a simulation, this function saves all parameter choices
save_setting <- function(num.trees, beta,
                         B, predictvar, sim, save_filename,
                         subdir = "Results",
                         reps, x,
                         num.features,
                         honesty,
                         response.scaling,
                         compute.oob.predictions,
                         min.node.size,
                         sampling,
                         data_foo,
                         n_vec,
                         alpha_ks,
                         d, p, q, get_truth, mu_diff_truth = NULL
) {
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

  setup <- list(num.trees = num.trees, beta = beta,
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
                d = d, p = p, q = q,
                get_truth = get_truth,
                mu_diff_truth = mu_diff_truth)
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
data_foo_cate <- function(sim, n, q, p, d, x = NULL) {
  #' p = dimension of X
  #' d = dimension of Y
  #' q = dimension of theta

  X <- if (is.null(x)) {
    matrix(runif(p * n), ncol = p, nrow = n, byrow = T)
  } else {
    rep(x, n) %>%
      matrix(nrow = n, byrow = TRUE)
  }

  if (sim == 1) { # First simulation (confounder affecting W and Y)
    stopifnot(q == 1)
    stopifnot(d == 2)
    stopifnot(p >= 3)
    W <- rbinom(n, 1, (1 + dbeta(X[, 3], 2, 4)) / 4)
    Y <-
      cbind(2 * (X[, 3] - 0.5) + rnorm(n))
    y <- seq(from = -3.5, to = 4, length.out = 50)

  } else if (sim == 2) { # Second simulation (heterogenous treatment effect)
    stopifnot(d == 2)
    stopifnot((q == 1) | (q == 2))
    W <- rbinom(n, 1, 0.5)
    etaX <- apply(X, 2, function(x) fooY(x))
    Y <-
      cbind((W - 1 / 2) * etaX[, 1] * etaX[, 2] + rnorm(n))
    y <- seq(from = -5.5, to = 5, length.out = 50)

  } else if (sim == 3) { # Third simulation (heterogeneous treatment effect + confounding)
    stopifnot(d == 2)
    stopifnot((q == 1) | (q == 2))
    stopifnot(p >= 3)
    Y <- matrix(rnorm(d * n), ncol = d, nrow = n, byrow = T)
    px3 <- 1 / 4 * (1 + dbeta(X[, 3], shape1 = 2, shape2 = 4))
    etaX <- apply(X, 2, function(x) fooY(x)) %>%
      rbind()
    # Y[, 2] is the treatment assignment W:
    W <- rbinom(n = n, size = 1, prob = px3)
    Y <- 2 * (X[, 3] - 1 / 2) - 1 / 2 * etaX[, 1] * etaX[, 2] +
      W * etaX[, 1] * etaX[, 2] + rnorm(n) %>%
      cbind()
    y <- seq(from = -5, to = 5, length.out = 50)

  } else if (sim == 4) {
    stopifnot(d == 2)
    stopifnot(p >= 2)
    stopifnot((q == 1) | (q == 2))
    W <- rbinom(n = n, size = 1, prob = expit(4 * X[, 2] - 2))
    Y <-
      cbind(100 * X[, 2] ^ 2 + (W - 0.5) * sin(3 * X[, 1]) + rnorm(n))
    y <- seq(from = -3, to = 90, length.out = 50)
  }

  list(X = X, Y = Y, W = W,
       x = get_truth_cate(sim = sim, q = q, p = p)$x,
       y = y)
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

# New data simulation scenarios for CoDiTE
data_foo_codite <- function(sim, n, q = 1, p = 1, d = 1, x = NULL) {
  X <- if (is.null(x)) {
    rnorm(n) * 0.5
  } else {
    rep(x, n)
  }
  if (is.null(x)) {
    x <- 0.5 %>% matrix()
  }

  W <- rbinom(n, size = 1, prob = exp(-X) / (1 + exp(-X)))

  if (sim == 1) {
    Y <- 2 * X + rnorm(n, mean = 1, sd = 0.5)
    y <- seq(from = -3, to = 5, length.out = 50) %>% cbind()
  } else if (sim == 2) {
    Y <- ifelse(W == 1, rnorm(n, mean = -1, sd = 1), rnorm(n, mean = 1, sd = 0.5))
    y <- seq(from = -5, to = 3, length.out = 50) %>% cbind()
  }
  return(list(X = matrix(X), W = matrix(W), Y = matrix(Y), x = x, y = y))
}

get_truth_codite_foo <- function(y, x, n_approx, sim, p, q, d, data_foo) {
  sample_data <- data_foo(sim = sim, n = n_approx, q = q, p = p, d = d, x = x)

  bandwidth_Y <- drf:::medianHeuristic(sample_data$Y)
  k_Y <- rbfdot(sigma = bandwidth_Y)

  K1_mean <-
    kernelMatrix(k_Y, sample_data$Y[sample_data$W == 1], y = y) %>%
    colMeans()
  K0_mean <-
    kernelMatrix(k_Y, sample_data$Y[sample_data$W == 0], y = y) %>%
    colMeans()
  K1_mean - K0_mean
}