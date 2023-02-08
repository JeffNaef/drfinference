# compute per-halfsample DRF weights and overall weights
bootdist <- function(X, Y, x, num.features, honesty, num.trees, B, beta,
                     response.scaling, compute.oob.predictions, min.node.size,
                     sampling = "binomial") {
  
  n <- dim(X)[1]
  nrowx <- nrow(x)
  
  # compute point estimator and DRF per halfsample (B halfsamples)
  # weightsb: B times n matrix of weights
  weightsb <- lapply(seq_len(B), function(b) {
    
    # half-sample index: according to theory, every datapoint is in or out
    # according to Bernoulli(0.5). In practice, may just take half of the
    # data (the else-option)
    indexb <- if (sampling == "binomial") {
      seq_len(n)[as.logical(rbinom(n, size = 1, prob = 0.5))]
    } else {
      sample(seq_len(n), floor(n / 2), replace = FALSE)
    }
    
    # fit DRF to half-sample
    DRFb <- drf(X = X[indexb, , drop = F],
                Y = Y[indexb, , drop = F],
                num.features = num.features,
                honesty = honesty,
                num.trees = num.trees,
                sample.fraction = (n) ^ beta / length(indexb), 
                response.scaling = response.scaling,
                compute.oob.predictions = compute.oob.predictions,
                min.node.size = min.node.size,
                ci.group.size = 1) 
    
    # extract weights (nrowx times number of split observations)
    weightsbfinal <- Matrix(0, nrow = nrowx, ncol = n, sparse = TRUE)
    weightsbfinal[, indexb] <- predict(DRFb, x)$weights 
    
    weightsbfinal
  })
  
  # postprocess weights. Every row of x gets its own list entry in weightsb.
  # Every list entry of weightsb contains a (B x n) matrix containing the
  # weights for all datapoints per halfsample DRF, and there are B many such
  # halfsample DRF's
  weightsb <- lapply(seq_len(nrowx), function(i) {
    do.call(rbind, sapply(seq_len(B), function(bb) {
      weightsb[[bb]][i, , drop = FALSE]
    }))
  })
  
  # compute overall weights as average from individual hslfsample DRF's
  weights_all <- lapply(seq_len(nrowx), function(i) {
    colMeans(weightsb[[i]]) %>%
      Matrix(sparse = T, nrow = 1)
  })
  
  # return
  list(weights_all = weights_all, weightsdist = weightsb)
}

# compute point estimator and variance
bootvar <- function(X, Y, param_foo, x,
                    num.features, honesty, num.trees, B, beta,
                    response.scaling, compute.oob.predictions, min.node.size,
                    sampling = "binomial") {
  nrowx <- nrow(x)
  
  # get per-halfsample (B many in total) and overall weights
  weights <-
    bootdist(X = X, Y = Y, x = x,
             num.features = num.features, honesty = honesty,
             num.trees = num.trees, B = B, beta = beta,
             response.scaling = response.scaling,
             compute.oob.predictions = compute.oob.predictions,
             min.node.size = min.node.size,
             sampling = sampling)
  
  # compute point estimator per halfsample.
  # every entry in the list corresponds to one row of x.
  thetahatb <- lapply(seq_len(nrowx), function(i) {
    sapply(seq_len(B), function(b) {
      param_foo(weights = weights$weightsdist[[i]][b, , drop = F],
                X = NULL, Y = Y, x = x[i, 1], probs = probs)
    }) %>%
      t()
  })
  
  # compute theta from the full DRF weights.
  # every entry in the list corresponds to one row of x.
  thetahat <-
    lapply(seq_len(nrowx), function(i) {
      param_foo(weights = weights$weights_all[[i]], X = X, Y = Y, x = x[i, ], probs = probs)
    })
  
  # compute variance from the halfsampling point estimators.
  # every entry in the list corresponds to one row of x.
  var_est <- if (length(thetahat[[1]]) > 1) {
    lapply(seq_len(nrowx), function(i) {
      a <- sweep(thetahatb[[i]], 2, thetahat[[i]], FUN = "-")
      crossprod(a, a) / B
    })
  } else {
    lapply(seq_len(nrowx), function(i) {
      mean((c(thetahatb[[i]]) - c(thetahat[[i]])) ^ 2)
    })
  } # end var_est
  
  # return
  list(thetahat = thetahat, var_est = var_est)
}
