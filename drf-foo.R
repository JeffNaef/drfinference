drfCI <- function(X, Y, B, sampling = "binomial", ...) {

  n <- dim(X)[1]

  # compute point estimator and DRF per halfsample
  # weightsb: B times n matrix of weights
  DRFlist <- lapply(seq_len(B), function(b) {

    # half-sample index
    indexb <- if (sampling == "binomial") {
      seq_len(n)[as.logical(rbinom(n, size = 1, prob = 0.5))]
    } else {
      sample(seq_len(n), floor(n / 2), replace = FALSE)
    }

    ## Using normal Bootstrap on the data and refitting DRF
    DRFb <-
      drf(X = X[indexb, , drop = F], Y = Y[indexb, , drop = F],
          ci.group.size = 1, ...)


    return(list(DRF = DRFb, indices = indexb))
  })
  
  return(list(DRFlist = DRFlist, X = X, Y = Y))
}


predictdrf <- function(DRF, x, functional = NULL, ...) {
  
  ntest <- nrow(x)
  n <- nrow(DRF$Y)
  
  
  weightsb <- lapply(DRF$DRFlist, function(l) {
    
    weightsbfinal <- Matrix(0, nrow = ntest, ncol = n, sparse = TRUE)
    
    weightsbfinal[, l$indices] <- predict(l$DRF, x)$weights 
    
    return(weightsbfinal)
  })
  
  
  weightsall <- Reduce("+", weightsb) / length(weightsb)
  
  if (!is.null(functional)) {
    stopifnot("Not yet implemented for several x" = ntest == 1)
    
    thetahatb <- 
      lapply(weightsb, function(w)  
        functional(weights = w, X = DRF$X, Y = DRF$Y, x = x))
    thetahatbforvar <- do.call(rbind, thetahatb)
    thetahat <- functional(weights = weightsall, X = DRF$X, Y = DRF$Y, x = x)
    thetahat <- matrix(thetahat, nrow = dim(x)[1])
    var_est <- if (dim(thetahat)[2] > 1) {  
      a <- sweep(thetahatbforvar, 2, thetahat, FUN = "-")
      crossprod(a, a) / B
    } else {
      mean((c(thetahatbforvar) - c(thetahat)) ^ 2) 
    }
    
    return(list(weights = weightsall, thetahat = thetahat, weightsb = weightsb, 
                var_est = var_est))
    
  } else {
    return(list(weights = weightsall, weightsb = weightsb))
  }
}
