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

do_plots <- function(n_vec = NULL, do_hist = TRUE, runout = 50) {
  alpha <- 0.05
  load("preliminary.RData")
  n_vec <- if (is.null(n_vec)) {
    setup$n_vec
  } else {
    n_vec
  }
  reps <- setup$reps
  z <- pnorm(1 - alpha / 2)

  save_filename <- setup$save_filename
  myColors <- 
    c("#117DE1", "#FC8D62", "#A6D854", "#FFDB6D", "#FF5733", "red", "cyan")[seq_len(length(n_vec))]
  names(myColors) <- n_vec
  
  # concatenate results for individual n's
  res_all <- tibble()
  for (n in n_vec) {
    print(paste0("n = ", n))
    
    res_new <- tibble()
    for (rr in seq_len(runout)) {
      if (runout == 0) {
        load(paste0(save_filename, "_n_", n, ".RData"))
      } else {
        load(paste0(save_filename, "_n_", n, "_run_", rr, ".RData"))
      }
      
      res_new_new <- res %>%
        mutate(n = n, 
               pvalue_small = (pvalue <= alpha)) 
      res_new <- rbind(res_new, res_new_new)
    }
    
    pvalue_alpha <- mean(res_new$pvalue_small)
    CIs <- get.CIs(pvalue_alpha, reps * runout, z)
    res_new <- res_new %>%
      mutate(pvalue_alpha = pvalue_alpha, 
             CI.lower = CIs[1], 
             CI.upper = CIs[2])
    res_all <- rbind(res_all, res_new)
  } # end for n
  
  # set line type
  my_alpha <- 0.15
  my_cex <- 2
  my_lwd <- 0.75 
  my_lwd_small <- 0.5
  textsize <- 10
  hadjust <- 0.5
  no_lines <- 1
  
  res_summary <- res_all %>%
    group_by(n) %>%
    slice(1)
  
  print(res_summary)
  dput(res_summary, file = "res_summary.txt")
  dput(res_summary, file = "res_summary.RData")
  
  p_alpha <- 
    res_summary %>%
    group_by(n) %>%
    ggplot(aes(x = n)) +
    geom_ribbon(aes(ymin = CI.lower, ymax = CI.upper), alpha = my_alpha, fill = "black") + 
    geom_hline(aes(yintercept = alpha), lwd = my_lwd_small, lty = 2) +
    geom_line(aes(y = pvalue_alpha), lwd = my_lwd, col = "black") + 
    geom_jitter(aes(y = pvalue_alpha), cex = my_cex, width = 0.00, height = 0.00, col = "black") + 
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
    labs(x = quote(n), title = "Proportion of\n p-values below 0.05", y = NULL)
  
  p_hist <- 
    res_all %>%
    ggplot(aes(x = pvalue)) + 
    geom_histogram(aes(y=..density..), colour = "gray", fill = "gray", bins = 30) + 
    #facet_wrap(~ n) + 
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
    labs(x = "observed p-values", title = "Histogram")
  
  p_qq <- 
    res_all %>%
    group_by(n) %>%
    ggplot(aes(sample = pvalue)) + 
    geom_qq(distribution = stats::qunif, cex = 0.5, color = "black") + 
    geom_qq_line(distribution = stats::qunif, color = "gray", lwd = 0.5) + 
    #facet_wrap(~ n) + 
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
    labs(x = "theoretical Uniform(0, 1)", title = "QQ-plot")
  
  p_dominated <- 
    res_all %>%
    group_by(n) %>%
    ggplot(aes(x = pvalue)) + 
    stat_ecdf(col = "black", cex = 0.9) + 
    geom_line(aes(y = pvalue), color = "gray", lwd = 0.5) + 
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
    labs(x = expression(p), title = "Empirical distribution function",
      y = expression(hat(F)(p)))
  
  
  num.plots <- length(n_vec)
  widths <- c(max((num.plots - 1), 1) / num.plots, 1 / num.plots)
  p_final <-
    ggarrange(p_qq, 
              p_hist,
              p_dominated, 
              ncol = 3, nrow = 1, common.legend = TRUE,
              legend = "right", widths = widths)
  annotate_figure(p_final, fig.lab = "Observed p-values", fig.lab.face = "bold", fig.lab.size = 14,
                  top = text_grob("")) 
  p_final
  
  ggsave("pvalues.pdf", width = 9, height = 3)
}

do_plots_CIbands <- function(n_vec = NULL, load_truth = FALSE, runout = 50) {
  load("preliminary.RData")
  n_vec <- if (is.null(n_vec)) {
    setup$n_vec
  } else {
    n_vec
  }
  
  reps <- setup$reps
  q <- setup$q
  alpha <- setup$alpha_ks
  z <- pnorm(1 - alpha / 2)
  mu_diff_truth <- if (load_truth) {
    load("mu_diff_truth.RData")
    mu_diff_truth
  } else {
    setup$mu_diff_truth
  }

  save_filename <- setup$save_filename
  myColors <- 
    c("#117DE1", "#FC8D62", "#A6D854", "#FFDB6D", "#FF5733", "red", "cyan")[seq_len(length(n_vec))]
  names(myColors) <- n_vec
  
  # concatenate results for individual n's
  res_all <- tibble()
  for (n in n_vec) {
    print(paste0("n = ", n))
    
    res_new <- tibble()
    for (rr in seq_len(runout)) {
      rept <- if (runout == 0) {
        load(paste0(save_filename, "_n_", n, ".RData"))
        rep(seq_len(reps), each = length(mu_diff_truth))
      } else {
        load(paste0(save_filename, "_n_", n, "_run_", rr, ".RData"))
        rep(seq_len(reps) + (rr - 1) * reps, each = length(mu_diff_truth))
      }
      
      res_new_new <- res %>%
        mutate(repetition = rept)
      res_new <- rbind(res_new, res_new_new)
    }
    
    cover <- mean(res_new$is.in)
    CIs <- get.CIs(cover, reps * runout, z)
    res_new <- res_new %>%
      mutate(cover = cover,
             CI.lower = CIs[1], 
             CI.upper = CIs[2])
    res_all <- rbind(res_all, res_new)
    
    
  } # end for n
  
  # set line type
  my_alpha <- 0.15
  my_cex <- 2
  my_lwd <- 0.75 
  my_lwd_small <- 0.5
  textsize <- 10
  hadjust <- 0.5
  no_lines <- 1
  
  res_summary <- res_all %>%
    group_by(n) %>%
    slice(1)
  
  print(res_summary %>% as.data.frame())
  dput(res_summary, file = "res_summary.txt")
  dput(res_summary, file = "res_summary.RData")
  
  p_cov <- res_summary %>%
    ggplot(aes(x = n)) +
    geom_ribbon(aes(ymin = CI.lower, ymax = CI.upper), alpha = my_alpha) + 
    geom_hline(aes(yintercept = 1 - alpha), lwd = my_lwd_small, lty = 2) +
    geom_line(aes(y = cover), lwd = my_lwd) + 
    geom_jitter(aes(y = cover), cex = my_cex, width = 0.00, height = 0.00) + 
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
         x = quote(n), title = "Simultaneous CI bands", y = NULL) #+
  #facet_wrap(~ name, ncol = 1)
  
  p_cov
  ggsave("confbands_coverage.pdf", width = 6, height = 5)
  
  true_func_dat <- 
    tibble(witness_truth = rep(mu_diff_truth, length(n_vec)), 
           n = rep(n_vec, each = length(mu_diff_truth)), 
           y = res_all %>% filter(repetition == 1) %>% pull(y))
  
  p_show_all <- res_all %>%
    group_by(n, repetition) %>%
    ggplot(aes(x = y)) + 
    geom_line(aes(y = CBl, group = repetition), color = "gray") + 
    geom_line(aes(y = CBu, group = repetition), color = "gray") + 
    geom_line(aes(y = hatmun_y, group = repetition), color = "#ADD8E6") + 
    geom_line(data = true_func_dat, aes(y = witness_truth), color = "black", lwd = 0.5) + 
    #facet_wrap(~ n) + 
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
         x = expression(y), title = "Coverage from CI bands", y = NULL) 
  
  p_show_all
  ggsave("confbands_plots.pdf", width = 4, height = 3)
}
