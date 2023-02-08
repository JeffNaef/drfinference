# compute confidence interval for estimated probability parameter (p.hat)
# of Bernoulli distribution that has been computed with rep.times
# many repetitions.
get.CIs <- function(p.hat, rep.times, z) {
  if (p.hat == 0) {
    ci.est <- c(0, 3 / rep.times)
  } else if (p.hat == 1) {
    ci.est <- c(1 - 3 / rep.times, 1)
  } else {
    ci.est <-
      c(p.hat - z * sqrt(p.hat * (1 - p.hat) / rep.times),
        p.hat + z * sqrt(p.hat * (1 - p.hat) / rep.times))
  }
  ci.est[1] <- pmax(ci.est[1], 0)
  ci.est[2] <- pmin(ci.est[2], 1)
  ci.est
}

# https://rpubs.com/huanfaChen/squash_remove_y_axix_ggplot
library(scales)
squash_axis <- function(from, to, factor) { 
  # A transformation function that squashes the range of [from, to] by factor on a given axis 
  
  # Args:
  #   from: left end of the axis
  #   to: right end of the axis
  #   factor: the compression factor of the range [from, to]
  #
  # Returns:
  #   A transformation called "squash_axis", which is capsulated by trans_new() function
  
  trans <- function(x) {
    # get indices for the relevant regions
    isq <- (x >= from) & (x <= to)
    ito <- (x > to)
    
    # apply transformation
    if (!any(is.na(isq))) {
      x[isq] <- from + (x[isq] - from) / factor
    }
    if (!any(is.na(ito))) {
      x[ito] <- from + (to - from) / factor + (x[ito] - to)
    }
    
    return(x)
  }
  
  inv <- function(x) {
    # get indices for the relevant regions
    isq <- (x >= from) & (x <= (from + (to - from) / factor))
    ito <- x > from + (to - from) / factor
    
    # apply transformation
    if (!any(is.na(isq))) {
      x[isq] <- from + (x[isq] - from) * factor
    }
    if (!any(is.na(ito))) {
      x[ito] <- to + (x[ito] - (from + (to - from) / factor))
    }
    
    return(x)
  }
  
  # return the transformation
  return(trans_new("squash_axis", trans, inv))
}

do_plots <- function(alpha = 0.05, n_vec = NULL, do_copula = FALSE, 
                     from = 0, to = 0.8, xoneproblematic_ind = c(11), factor = 20, 
                     from2 = 0, to2 = 0.8, factor2 = 20) {
  load("preliminary.RData")
  n_vec <- if (is.null(n_vec)) {
    setup$n_vec
  } else {
    n_vec
  }
  num.trees <- setup$num.trees
  beta <- setup$beta
  thetazero <- setup$thetazero
  reps <- setup$reps
  B <- setup$B
  q <- setup$q
  qch <- qchisq(1 - alpha, df = q)
  do_quantile <- setup$do_quantile
  do_copula <- setup$do_copula
  probs <- setup$probs
  save_filename <- setup$save_filename
  theta_names <- if (!do_quantile) {
    sapply(seq_len(q), function(qq) paste0("th", qq))
  } else {
    probs
  }
  x <- setup$x
  numx <- nrowx <- nrow(x)
  all_cols <- 
    c("#117DE1", "#FC8D62", "#A6D854", "#FFDB6D", "#FF5733", "red", "cyan")
  # set colors
  if (!do_quantile & !do_copula) {
    myColors <- all_cols[1:(do_CF + 1)] 
    names(myColors) <- c("DRF", "GRF")[c(TRUE, do_CF)]
  } else if (do_quantile) {
    myColors <- all_cols[seq_len(length(theta_names))] 
    names(myColors) <- theta_names
  } else if (do_copula) {
    myColors <- all_cols[seq_len(length(n_vec))] 
    names(myColors) <- as.character(n_vec)
  }
  
  # concatenate results for individual n's
  res_all <- tibble()
  for (n in n_vec) {
    print(paste0("n = ", n))
    load(paste0(save_filename, "_n_", n, ".RData"))
    
    ### DRF
    # point estimator
    thetahat <- do.call(rbind, res[, "thetahat"])
    thetahat <- lapply(seq_len(numx), function(i) {
      a <- do.call(rbind, thetahat[, i])
      rownames(a) <- NULL
      a
    })
    
    # variance
    a <- do.call(rbind, res[, "var_est"])
    var_est <- lapply(seq_len(numx), function(i) {
      aa <- if (q == 1) {
        a_interm <- do.call(rbind, a[, i])
        rownames(a_interm) <- NULL
        a_interm
      } else {
        sapply(seq_len(reps), function(rr) {
          diag(a[, i][[rr]])
        }) %>%
          t()
      } # end var_est
      aa
    })
    
    # confidence intervals
    qnormalph <- qnorm(1 - alpha / 2)
    CI.lower <- lapply(seq_len(numx), function(j) {
      aa <- sapply(seq_len(q), function(i) {
        th <- thetahat[[j]][, i]
        sd_th <- sqrt(var_est[[j]][, i])
        th - qnormalph * sd_th
      })
      aa
    })
    CI.upper <- lapply(seq_len(numx), function(j) {
      aa <- sapply(seq_len(q), function(i) {
        th <- thetahat[[j]][, i]
        sd_th <- sqrt(var_est[[j]][, i])
        th + qnormalph * sd_th
      })
      aa
    })
    
    # ellipse CI
    # Z~N(0,I), P(Z \in \{x: x'x <= q \} ) = P(Z'Z <= q) = P(chi^2_d <=q), so we need qch=qchisq(0.95, df=d)
    is.in.ell <- mclapply(seq_len(numx), function(j) {
      thz <- if (!do_quantile & !do_copula) {
        thetazero
      } else {
        thetazero[j, ]
      }
      aa <- sapply(seq_len(reps), function(i) {
        th <- thetahat[[j]][i, ]
        sd_mat <- solve(sqrtm(as.matrix(a[, j][[i]])))
        th_scaled <- (th - thz) %*% sd_mat
        tcrossprod(th_scaled) <= qch
      })
      aa
    }, mc.cores = 12)
    
    # concatenate results
    res_all_new <-
      tibble(thhat = do.call(c, thetahat),
             var = do.call(c, var_est),
             name = rep(rep(theta_names, each = reps), numx),
             n = n,
             thzero = if (!do_quantile & !do_copula) {
               rep(rep(thetazero, each = reps), numx)
             } else {
               do.call(c, lapply(seq_len(nrowx), function(i) sapply(seq_len(reps), function(j) thetazero[i, ]) %>% t()))},
             ci.lower = do.call(c, CI.lower),
             ci.upper = do.call(c, CI.upper),
             is.in = (ci.lower <= thzero) & (thzero <= ci.upper),
             bias = thhat - thzero,
             ci.len = ci.upper - ci.lower,
             is.in.ellipse = do.call(c, lapply(seq_len(numx), function(i) {
               rep(is.in.ell[[i]], q)
             })),
             method = "DRF",
             xone = rep(x[, 1], each = q * reps))
    
    res_all <- rbind(res_all, res_all_new)
    
    # incorporate grf results (only works for one-dimensional theta)
    if (do_CF) {
      stopifnot(nrowx == 1)
      thetahat_CF <- do.call(rbind, res[, "thetahat_CF"])
      rownames(thetahat_CF) <- NULL
      
      var_est_CF <- do.call(rbind, res[, "var_est_CF"])
      rownames(var_est_CF) <- NULL
      
      CI.lower_CF <- sapply(seq_len(reps), function(i) {
        th <- thetahat_CF[i, ]
        sd_th <- sqrt(var_est_CF[i, ])
        th - qnorm(1 - alpha / 2) * sd_th
      })
      CI.upper_CF <- sapply(seq_len(reps), function(i) {
        th <- thetahat_CF[i, ]
        sd_th <- sqrt(var_est_CF[i, ])
        th + qnorm(1 - alpha / 2) * sd_th
      })
      res_new_CF <-
        tibble(thhat = c(thetahat_CF),
               var = c(var_est_CF),
               name = rep(theta_names, each = reps),
               n = n,
               thzero = rep(thetazero, each = reps),
               ci.lower = CI.lower_CF,
               ci.upper = CI.upper_CF,
               is.in = (ci.lower <= thzero) & (thzero <= ci.upper),
               bias = thhat - thzero,
               ci.len = ci.upper - ci.lower,
               is.in.ellipse = FALSE,
               method = "GRF",
               xone = rep(x[1, 1], reps * q))
      res_all <- rbind(res_all, res_new_CF)
    }
    
  }
  
  xone_all <- res_all %>%
    distinct(xone) %>%
    pull(xone) %>%
    sort()
  
  res_all <- res_all %>%
    group_by(method, n, name, xone) %>%
    mutate(cover = mean(is.in),
           cover.ell = mean(is.in.ellipse),
           sd_overall = mean(sqrt(var)),
           th_mean_overall = mean(thhat),
           scale_thhat = (thhat - thzero) / sd_overall,
           sd_scale_thhat = sd(scale_thhat),
           CI_len_med = median(ci.len),
           bias_med = median(bias),
           sd_true = sd(thhat))
  
  
  ##################### plot a histogram
  
  mean_val_dat <- res_all %>%
    group_by(method, n, name) %>%
    dplyr::select(method, n, name, scale_thhat) %>%
    mutate(mean_val = mean(scale_thhat),
           sd_val = sd(scale_thhat))
  
  approx_tib <-
    tibble(xx = c(res_all$scale_thhat, res_all$scale_thhat),
           gauss_approx = c(dnorm(res_all$scale_thhat,
                                  mean = mean_val_dat$mean_val,
                                  sd = mean_val_dat$sd_val),
                            dnorm(res_all$scale_thhat, mean = 0, sd = 1)),
           type = c(rep("empirical var", nrow(res_all)), rep("N(0, 1)", nrow(res_all))),
           name = c(res_all$name, res_all$name),
           n = c(res_all$n, res_all$n),
           method = c(res_all$method, res_all$method),
           mean_val = c(mean_val_dat$mean_val, rep(0, nrow(res_all)))) %>%
    group_by(method, n, name)
  
  if (nrowx <= 3) {
    res_all %>%
      ggplot() +
      geom_histogram(aes(y = ..density.., x = scale_thhat), bins = 80) +
      geom_line(data = approx_tib, aes(x = xx, y = gauss_approx, col = type, lty = type),
                lwd = 0.75) +
      geom_vline(data = approx_tib, aes(xintercept = mean_val, col = type, lty = type),
                 show.legend = FALSE, lwd = 0.75) +
      ggtitle(paste0("B = ", B, ", reps = ", reps, ", trees = ", num.trees)) +
      xlab("(theta^ - theta0) / sd^(theta^)") +
      ylab("") +
      facet_wrap(method + name ~ n + xone)
    
    ggsave(filename = "hist.pdf")
  }
  
  
  ##################### coverage plot
  
  # set line type
  my_alpha <- 0.15
  my_cex <- 0.7
  my_lwd <- 0.3 
  my_lwd_small <- 0.3
  textsize <- 10
  hadjust <- 0.5
  no_lines <- 1
  
  # confidence bands for coverage plot
  CI_cov <- sapply(seq_len(nrow(res_all)), function(i) {
    get.CIs(res_all$cover[i], reps, pnorm(1 - alpha / 2))
  }) %>%
    t()
  CI_cov_ell <- sapply(seq_len(nrow(res_all)), function(i) {
    get.CIs(res_all$cover.ell[i], reps, pnorm(1 - alpha / 2))
  }) %>%
    t()
  res_all <- res_all %>%
    ungroup() %>%
    mutate(CI_cov_low = CI_cov[, 1],
           CI_cov_upp = CI_cov[, 2],
           CI_cov_low_ell = CI_cov_ell[, 1],
           CI_cov_upp_ell = CI_cov_ell[, 2])
  
  dat_cov <- res_all %>%
    group_by(method, n, name, xone) %>%
    slice(1) %>%
    mutate(name = as.character(name))
  if (do_copula) {
    dat_cov <- dat_cov %>%
      mutate(name = factor(n, levels = sort(n_vec)))
  }
  name_color <- ifelse(do_copula, "n", "Quantile")
  p_cov <- if (!do_quantile & !do_copula) {
    dat_cov %>%
      group_by(method, n, name, xone) %>%
      ggplot(aes(x = n)) +
      geom_ribbon(aes(ymin = CI_cov_low, ymax = CI_cov_upp, fill = method), alpha = my_alpha) + 
      geom_hline(aes(yintercept = 1 - alpha), lwd = my_lwd_small, lty = 2) +
      geom_line(aes(y = cover, col = method), lwd = my_lwd) + 
      geom_jitter(aes(y = cover, col = method), cex = my_cex, width = 0.00, height = 0.00) + 
      scale_colour_manual(name = "Method", values = myColors) +
      scale_fill_manual(name = "Method", values = myColors) +
      guides(colour = guide_legend(override.aes = list(alpha = 1))) +
      theme(plot.title = element_text(face = "bold", size = textsize, hjust = hadjust),
            legend.background = element_rect(fill = NA, size = 4, colour = NA),
            legend.key = element_blank(),
            plot.background = element_rect(fill = "white", colour = "white"),
            panel.background = element_rect(fill = "white", colour = "white"),
            axis.ticks = element_line(colour = "grey70", size = 0.2),
            panel.grid.major = element_line(colour = "grey70", size = 0.2),
            panel.grid.minor = element_line(colour = "grey70", size = 0.05),
            axis.text.x = element_text(angle = 90),
            panel.spacing = unit(no_lines, "lines"),
            strip.text.y = element_blank()) +
      labs(color = "Method", fill = "Method",
           x = quote(n), title = "Coverage", y = NULL) #+
    #facet_wrap(~ name, ncol = 1)
  } else if (do_quantile) {
    xoneproblematic <- xone_all[xoneproblematic_ind]
    bb <- c(seq(from, to, 0.2), seq(to, 1, 0.05))
    dat_low <- dat_cov %>%
      dplyr::filter(xone < min(xoneproblematic))
    dat_problematic <- dat_cov %>%
      dplyr::filter(xone %in% xoneproblematic)
    dat_high <- dat_cov %>%
      dplyr::filter(xone > max(xoneproblematic))
    
    dat_cov %>%
      group_by(method, n, name, xone) %>%
      ggplot(aes(x = xone)) +
      geom_hline(aes(yintercept = 1 - alpha), lwd = my_lwd_small, lty = 2) +
      geom_jitter(aes(y = cover, col = name), cex = my_cex, width = 0.00, height = 0.00) + 
      scale_colour_manual(name = name_color, values = myColors) +
      scale_fill_manual(name = name_color, values = myColors) +
      guides(colour = guide_legend(override.aes = list(alpha = 1))) +
      scale_y_continuous(trans = squash_axis(from, to, factor), 
                         breaks = bb, labels = bb,
                         limits = c(0, 1)) + 
      theme(plot.title = element_text(face = "bold", size = textsize, hjust = hadjust),
            legend.background = element_rect(fill = NA, size = 4, colour = NA),
            legend.key = element_blank(),
            plot.background = element_rect(fill = "white", colour = "white"),
            panel.background = element_rect(fill = "white", colour = "white"),
            axis.ticks = element_line(colour = "grey70", size = 0.2),
            panel.grid.major = element_line(colour = "grey70", size = 0.2),
            panel.grid.minor = element_line(colour = "grey70", size = 0.05),
            axis.text.x = element_text(angle = 90),
            panel.spacing = unit(no_lines, "lines"),
            strip.text.y = element_blank()) +
      labs(color = name_color, fill = name_color, 
           x = quote(x[1]), title = "Coverage", y = NULL) + 
      geom_line(data = dat_low, aes(y = cover, col = name), lwd = my_lwd) + 
      geom_line(data = dat_problematic, aes(y = cover, col = name), lwd = my_lwd) +
      geom_line(data = dat_high, aes(y = cover, col = name), lwd = my_lwd) + 
      geom_ribbon(data = dat_low, aes(ymin = CI_cov_low, ymax = CI_cov_upp, fill = name), alpha = my_alpha) + 
      geom_ribbon(data = dat_problematic, aes(ymin = CI_cov_low, ymax = CI_cov_upp, fill = name), alpha = my_alpha) + 
      geom_ribbon(data = dat_high, aes(ymin = CI_cov_low, ymax = CI_cov_upp, fill = name), alpha = my_alpha) #+
    #facet_wrap(~ n, ncol = 1)
  } else if (do_copula) {
    start_val <- 0.9
    start_val_low <- 0.1
    bb <- c(seq(start_val_low, start_val, 0.2), seq(start_val, 1, 0.05))
    dat_cov %>%
      dplyr::filter(cover >= start_val) %>% 
      group_by(method, n, name, xone) %>%
      ggplot(aes(x = xone)) +
      geom_ribbon(aes(ymin = CI_cov_low, ymax = CI_cov_upp, fill = name), alpha = my_alpha) + 
      geom_hline(aes(yintercept = 1 - alpha), lwd = my_lwd_small, lty = 2) +
      geom_line(aes(y = cover, col = name), lwd = my_lwd) + 
      geom_jitter(aes(y = cover, col = name), cex = my_cex, width = 0.00, height = 0.00) + 
      geom_point(data = dat_cov %>% dplyr::filter(cover < start_val), 
                 aes(x = xone, y = cover, color = name), cex = my_cex) + 
      scale_colour_manual(name = name_color, values = myColors) +
      scale_fill_manual(name = name_color, values = myColors) +
      guides(colour = guide_legend(override.aes = list(alpha = 1))) +
      theme(plot.title = element_text(face = "bold", size = textsize, hjust = hadjust),
            legend.background = element_rect(fill = NA, size = 4, colour = NA),
            legend.key = element_blank(),
            plot.background = element_rect(fill = "white", colour = "white"),
            panel.background = element_rect(fill = "white", colour = "white"),
            axis.ticks = element_line(colour = "grey70", size = 0.2),
            panel.grid.major = element_line(colour = "grey70", size = 0.2),
            panel.grid.minor = element_line(colour = "grey70", size = 0.05),
            axis.text.x = element_text(angle = 90),
            panel.spacing = unit(no_lines, "lines"),
            strip.text.y = element_blank()) +
      labs(color = name_color, fill = name_color, 
           x = quote(x[1]), title = "Coverage", y = NULL) + 
      scale_y_continuous(trans = squash_axis(0, start_val, 20), 
                         breaks = bb, labels = bb, 
                         limits = c(0, 1))
  }
  
  if (!do_quantile & !do_copula) {
    fac_rate <- dat_cov %>%
      filter(method == "DRF") %>%
      group_by(name) %>%
      dplyr::select(n, name, method, CI_len_med) %>%
      mutate(factor = mean(CI_len_med / n ^ - (1 - beta)))
    dat_convergence_rate <-
      tibble(n = rep(n_vec, q),
             CI_len_med = rep(n_vec ^ - (1 - beta) * fac_rate$factor * 0.9, q),
             method = "theoretic",
             name = rep(as.character(theta_names), each = length(n_vec)))
    dat_len <- bind_rows(dat_cov, dat_convergence_rate)
    
    p_sd <-
      dat_cov %>%
      ggplot(aes(x = n, lty = method)) +
      geom_point(aes(y = sd_overall), col = "red") +
      geom_line(aes(y = sd_overall), col = "red") +
      geom_point(aes(y = sd_true), col = "blue") +
      geom_line(aes(y = sd_true), col = "blue") +
      facet_wrap(~ name)
    p_sd
  }
  
  myCol_len <- myColors
  myCol_len <- c(myCol_len, "black")
  names(myCol_len) <- c(names(myColors), "theoretic")
  p_len <- if (!do_quantile & !do_copula) {
    dat_len %>%
      dplyr::filter(method != "theoretic") %>% 
      ggplot(aes(x = log10(n))) +
      scale_y_continuous(trans = "log10") +
      geom_line(aes(y = CI_len_med, col = method), lwd = my_lwd) + 
      geom_point(aes(y = CI_len_med, col = method), cex = my_cex) + 
      scale_color_manual(name = "Method", values = myCol_len) +
      scale_fill_manual(name = "Method", values = myCol_len) +
      guides(colour = guide_legend(override.aes = list(alpha = 1))) +
      theme(plot.title = element_text(face = "bold", size = textsize, hjust = hadjust),
            legend.background = element_rect(fill = NA, size = 4, colour = NA),
            legend.key = element_blank(),
            plot.background = element_rect(fill = "white", colour = "white"),
            panel.background = element_rect(fill = "white", colour = "white"),
            axis.ticks = element_line(colour = "grey70", size = 0.2),
            panel.grid.major = element_line(colour = "grey70", size = 0.2),
            panel.grid.minor = element_line(colour = "grey70", size = 0.05),
            axis.text.x = element_text(angle = 90),
            strip.text.y = element_blank(),
            panel.spacing = unit(no_lines, "lines")) +
      labs(color = "Method", fill = "Method", 
           x = quote(log(n)), title = "Log median CI length", y = NULL) #+
    #facet_wrap(~ name, ncol = 1)
  } else if (do_quantile) {
    dat_cov %>%
      ggplot(aes(x = xone)) +
      scale_y_continuous(trans = "log10") +
      geom_point(aes(y = CI_len_med, col = name), cex = my_cex) +
      scale_color_manual(name = name_color, values = myCol_len) +
      scale_fill_manual(name = name_color, values = myCol_len) +
      guides(colour = guide_legend(override.aes = list(alpha = 1))) +
      theme(plot.title = element_text(face = "bold", size = textsize, hjust = hadjust),
            legend.background = element_rect(fill = NA, size = 4, colour = NA),
            legend.key = element_blank(),
            plot.background = element_rect(fill = "white", colour = "white"),
            panel.background = element_rect(fill = "white", colour = "white"),
            axis.ticks = element_line(colour = "grey70", size = 0.2),
            panel.grid.major = element_line(colour = "grey70", size = 0.2),
            panel.grid.minor = element_line(colour = "grey70", size = 0.05),
            axis.text.x = element_text(angle = 90),
            strip.text.y = element_blank(),
            panel.spacing = unit(no_lines, "lines")) +
      labs(color = name_color, fill = name_color, 
           x = quote(x[1]), title = "Log median CI length", y = NULL) + 
      geom_line(data = dat_low, aes(y = CI_len_med, col = name), lwd = my_lwd) + 
      geom_line(data = dat_problematic, aes(y = CI_len_med, col = name), lwd = my_lwd) +
      geom_line(data = dat_high, aes(y = CI_len_med, col = name), lwd = my_lwd) 
  } else if (do_copula) {
    dat_cov %>%
      ggplot(aes(x = xone)) +
      scale_y_continuous(trans = "log10") +
      geom_line(aes(y = CI_len_med, col = name), lwd = my_lwd) +
      geom_point(aes(y = CI_len_med, col = name), cex = my_cex) + 
      scale_color_manual(name = name_color, values = myCol_len) +
      scale_fill_manual(name = name_color, values = myCol_len) +
      guides(colour = guide_legend(override.aes = list(alpha = 1))) +
      theme(plot.title = element_text(face = "bold", size = textsize, hjust = hadjust),
            legend.background = element_rect(fill = NA, size = 4, colour = NA),
            legend.key = element_blank(),
            plot.background = element_rect(fill = "white", colour = "white"),
            panel.background = element_rect(fill = "white", colour = "white"),
            axis.ticks = element_line(colour = "grey70", size = 0.2),
            panel.grid.major = element_line(colour = "grey70", size = 0.2),
            panel.grid.minor = element_line(colour = "grey70", size = 0.05),
            axis.text.x = element_text(angle = 90),
            strip.text.y = element_blank(),
            panel.spacing = unit(no_lines, "lines")) +
      labs(color = name_color, fill = name_color, 
           x = quote(x[1]), title = "Log median CI length", y = NULL)  
  }
 
  p_bias <- if (!do_quantile & !do_copula) {
    dat_cov %>%
      ggplot(aes(x = n)) +
      geom_hline(aes(yintercept = 0), lwd = my_lwd_small, lty = 2) +
      geom_line(aes(y = bias_med, col = method), lwd = my_lwd) + 
      geom_point(aes(y = bias_med, col = method), cex = my_cex) + 
      scale_color_manual(name = "Method", values = myColors) +
      scale_fill_manual(name = "Method", values = myColors) +
      guides(colour = guide_legend(override.aes = list(alpha = 1))) +
      theme(plot.title = element_text(face = "bold", size = textsize, hjust = hadjust),
            legend.background = element_rect(fill = NA, size = 4, colour = NA),
            legend.key = element_blank(),
            plot.background = element_rect(fill = "white", colour = "white"),
            panel.background = element_rect(fill = "white", colour = "white"),
            axis.ticks = element_line(colour = "grey70", size = 0.2),
            panel.grid.major = element_line(colour = "grey70", size = 0.2),
            panel.grid.minor = element_line(colour = "grey70", size = 0.05),
            axis.text.x = element_text(angle = 90),
            strip.text.y = element_blank(),
            panel.spacing = unit(no_lines, "lines")) +
      labs(color = "Method", fill = "Method", 
           x = quote(n), title = "Median bias", y = NULL) #+
    #facet_wrap(~ name, ncol = 1)
  } else if (do_quantile) {
    bb <- c(seq(-0.2, 0.2, 0.1), seq(0.2, 0.6, 0.2))
    dat_cov %>%
      ggplot(aes(x = xone)) +
      geom_hline(aes(yintercept = 0), lwd = my_lwd_small, lty = 2) +
      geom_point(aes(y = bias_med, col = name), cex = my_cex) + 
      scale_color_manual(name = name_color, values = myColors) +
      scale_fill_manual(name = name_color, values = myColors) +
      guides(colour = guide_legend(override.aes = list(alpha = 1))) +
      theme(plot.title = element_text(face = "bold", size = textsize, hjust = hadjust),
            legend.background = element_rect(fill = NA, size = 4, colour = NA),
            legend.key = element_blank(),
            plot.background = element_rect(fill = "white", colour = "white"),
            panel.background = element_rect(fill = "white", colour = "white"),
            axis.ticks = element_line(colour = "grey70", size = 0.2),
            panel.grid.major = element_line(colour = "grey70", size = 0.2),
            panel.grid.minor = element_line(colour = "grey70", size = 0.05),
            axis.text.x = element_text(angle = 90),
            strip.text.y = element_blank(),
            panel.spacing = unit(no_lines, "lines")) +
      labs(color = name_color, fill = name_color, 
           x = quote(x[1]), title = "Median bias", y = NULL) + 
      geom_line(data = dat_low, aes(y = bias_med, col = name), lwd = my_lwd) + 
      geom_line(data = dat_problematic, aes(y = bias_med, col = name), lwd = my_lwd) +
      geom_line(data = dat_high, aes(y = bias_med, col = name), lwd = my_lwd)  + 
      scale_y_continuous(trans = squash_axis(0.2, 0.6, 10), 
                         breaks = bb, labels = bb#,
      ) #+ 
    #facet_wrap(~ n, ncol = 1)
  } else if (do_copula) {
    dat_cov %>%
      ggplot(aes(x = xone)) +
      geom_hline(aes(yintercept = 0), lwd = my_lwd_small, lty = 2) +
      geom_line(aes(y = bias_med, col = name), lwd = my_lwd) +
      geom_point(aes(y = bias_med, col = name), cex = my_cex) + 
      scale_color_manual(name = name_color, values = myColors) +
      scale_fill_manual(name = name_color, values = myColors) +
      guides(colour = guide_legend(override.aes = list(alpha = 1))) +
      theme(plot.title = element_text(face = "bold", size = textsize, hjust = hadjust),
            legend.background = element_rect(fill = NA, size = 4, colour = NA),
            legend.key = element_blank(),
            plot.background = element_rect(fill = "white", colour = "white"),
            panel.background = element_rect(fill = "white", colour = "white"),
            axis.ticks = element_line(colour = "grey70", size = 0.2),
            panel.grid.major = element_line(colour = "grey70", size = 0.2),
            panel.grid.minor = element_line(colour = "grey70", size = 0.05),
            axis.text.x = element_text(angle = 90),
            strip.text.y = element_blank(),
            panel.spacing = unit(no_lines, "lines")) +
      labs(color = name_color, fill = name_color, 
           x = quote(x[1]), title = "Median bias", y = NULL) 
  }
  
  p_final <-
    ggarrange(p_cov, p_len, p_bias, ncol = 3, nrow = 1, common.legend = TRUE,
              legend = "right")
  p_final
  
  ggsave("coverage.pdf", width = 6, height = 3)
  
  # plot point estimators
  if (do_quantile | do_copula) {
    pp <- dat_cov %>%
      ggplot(aes(x = xone)) +
      geom_point(aes(y = thhat, col = name)) +
      geom_line(aes(y = thhat, col = name))
    pp <- if (do_quantile) {
      pp +
        facet_wrap(~ n, ncol = 1) +
        geom_line(aes(y = thzero, col = name))
    } else {
      pp +
        geom_line(aes(y = thzero), col = "black")
    }
    pp
    
    ggsave("estimated_params.pdf", width = 6, height = 4)
  }
  
  # plot ellipse information
  dat_ell <- if (!do_quantile & !do_copula) {
    res_all %>%
      filter(name == theta_names[1], method %in% c("DRF", "DRF-center")) %>%
      dplyr::select(n, name, method, cover.ell, CI_cov_low_ell, CI_cov_upp_ell, xone) %>%
      group_by(n, method) %>%
      slice(1) %>%
      ungroup()
  } else {
    res_all %>%
      filter(name == theta_names[1], method %in% c("DRF", "DRF-center")) %>%
      dplyr::select(n, name, method, cover.ell, CI_cov_low_ell, CI_cov_upp_ell, xone) %>%
      group_by(n, method, xone) %>%
      slice(1) %>%
      ungroup()
  }
  if (do_quantile) {
    bb <- c(seq(from2, to2, 0.2), seq(to2, 1, 0.05))
    dat_low <- dat_ell %>%
      filter(xone < xoneproblematic)
    dat_problematic <- dat_ell %>%
      filter(xone %in% xoneproblematic)
    dat_high <- dat_ell %>%
      filter(xone > xoneproblematic)
  }
  
  if (!do_quantile & !do_copula) {
    dat_ell %>%
      group_by(n) %>%
      ggplot(aes(x = n)) +
      geom_ribbon(aes(ymin = CI_cov_low_ell, ymax = CI_cov_upp_ell, fill = method), alpha = my_alpha) + #, pch = power_foo+
      geom_hline(aes(yintercept = 1 - alpha), lwd = my_lwd_small, lty = 2) +
      geom_line(aes(y = cover.ell, col = method), lwd = my_lwd) + # , lty = power_foo
      geom_jitter(aes(y = cover.ell, col = method), cex = my_cex, width = 0.00, height = 0.00) + #, pch = power_foo
      scale_colour_manual(name = "Method", values = myColors) +
      scale_fill_manual(name = "Method", values = myColors) +
      guides(colour = guide_legend(override.aes = list(alpha = 1))) +
      theme(plot.title = element_text(face = "bold", size = textsize, hjust = hadjust),
            legend.background = element_rect(fill = NA, size = 4, colour = NA),
            legend.key = element_blank(),
            plot.background = element_rect(fill = "white", colour = "white"),
            panel.background = element_rect(fill = "white", colour = "white"),
            axis.ticks = element_line(colour = "grey70", size = 0.2),
            panel.grid.major = element_line(colour = "grey70", size = 0.2),
            panel.grid.minor = element_line(colour = "grey70", size = 0.05),
            axis.text.x = element_text(angle = 90),
            panel.spacing = unit(no_lines, "lines"),
            strip.text.y = element_blank()) +
      labs(color = "Method", fill = "Method", 
           x = quote(n), title = "Coverage for ellipse", y = NULL)
  } else if (do_quantile) {
    
    bb <- c(seq(0.1, 0.9, 0.2), seq(0.9, 1, 0.05))
    
    dat_ell %>%
      group_by(n) %>%
      ggplot(aes(x = xone)) +
      geom_hline(aes(yintercept = 1 - alpha), lwd = my_lwd_small, lty = 2) +
      geom_jitter(aes(y = cover.ell), cex = my_cex, width = 0.00, height = 0.00) +
      scale_colour_manual(name = name_color, values = myColors) +
      scale_fill_manual(name = name_color, values = myColors) +
      guides(colour = guide_legend(override.aes = list(alpha = 1))) +
      theme(plot.title = element_text(face = "bold", size = textsize, hjust = hadjust),
            legend.background = element_rect(fill = NA, size = 4, colour = NA),
            legend.key = element_blank(),
            plot.background = element_rect(fill = "white", colour = "white"),
            panel.background = element_rect(fill = "white", colour = "white"),
            axis.ticks = element_line(colour = "grey70", size = 0.2),
            panel.grid.major = element_line(colour = "grey70", size = 0.2),
            panel.grid.minor = element_line(colour = "grey70", size = 0.05),
            axis.text.x = element_text(angle = 90),
            panel.spacing = unit(no_lines, "lines"),
            strip.text.y = element_blank()) +
      #facet_wrap(~ n, ncol = 1) +
      labs(color = "Method", fill = "Method", 
           x = quote(x[1]), title = "Coverage for ellipse CI's", y = NULL) + 
      scale_y_continuous(trans = squash_axis(from2, 0.9, 20), 
                         breaks = bb, labels = bb, 
                         limits = c(0.1, 1)) + 
      geom_line(data = dat_low, aes(y = cover.ell), lwd = my_lwd) + 
      geom_line(data = dat_problematic, aes(y = cover.ell), lwd = my_lwd) +
      geom_line(data = dat_high, aes(y = cover.ell), lwd = my_lwd) + 
      geom_ribbon(data = dat_low, aes(ymin = CI_cov_low_ell, ymax = CI_cov_upp_ell), alpha = my_alpha) + 
      geom_ribbon(data = dat_problematic, aes(ymin = CI_cov_low_ell, ymax = CI_cov_upp_ell), alpha = my_alpha) + 
      geom_ribbon(data = dat_high, aes(ymin = CI_cov_low_ell, ymax = CI_cov_upp_ell), alpha = my_alpha)
  }
  
  ggsave(filename = "coverage-ellipse.pdf", width = 5, height = 3)
}
