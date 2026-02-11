gen_logn <- function(mu, sigma, time_adj, p_val, n_simu, seed) {
  
  p_val    <- as.numeric(p_val)
  sigma    <- as.numeric(sigma)
  time_adj <- as.numeric(time_adj)
  
  width_p <- length(p_val)
  width_s <- length(sigma)
  width_t <- length(time_adj)
  width_v <- min(width_p, width_s, width_t)
  
  p_val    <- p_val[seq_len(width_v)]
  sigma    <- sigma[seq_len(width_v)]
  time_adj <- time_adj[seq_len(width_v)]
  
  sigmaT  <- sqrt(time_adj) * sigma
  mean_ZT <- (mu - (sigma^2) / 2) * time_adj
  
  ln_p <- log(p_val)
  
  ln_p_mat   <- matrix(rep(ln_p,    each = n_simu), nrow = n_simu)
  mean_mat   <- matrix(rep(mean_ZT, each = n_simu), nrow = n_simu)
  sigmaT_mat <- matrix(rep(sigmaT,  each = n_simu), nrow = n_simu)
  
  set.seed(seed)
  rv_norm <- matrix(rnorm(n_simu * width_v), nrow = n_simu)
  
  ln_pr  <- ln_p_mat + mean_mat + rv_norm * sigmaT_mat
  p_simu <- exp(ln_pr)
  
  p_dev <- (p_simu - exp(ln_p_mat)) / exp(ln_p_mat)
  
  diff_m <- p_val - colMeans(p_simu)
  p_simu <- p_simu + matrix(rep(diff_m, each = n_simu), nrow = n_simu)
  
  list(p_dev = p_dev, p_simu = p_simu)
}