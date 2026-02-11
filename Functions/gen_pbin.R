gen_pbin <- function(p_simu, bins, idx_bin) {
  
  p_simu <- as.matrix(p_simu)
  bins   <- as.matrix(bins)
  
  if (idx_bin == 1) {
    # one bin vector for all crops
    bin_s   <- bins[, 1, drop = TRUE]
    width_s <- ncol(p_simu)
    
    # replicate bin_s across columns (same as kron(ones(1,width_s),bin_s))
    bin_mat <- matrix(rep(bin_s, times = width_s),
                      nrow = length(bin_s), ncol = width_s, byrow = FALSE)
    
  } else {
    width_s <- min(ncol(p_simu), ncol(bins))
    bin_mat <- bins[, seq_len(width_s), drop = FALSE]
    p_simu  <- p_simu[, seq_len(width_s), drop = FALSE]
  }
  
  p_bins <- matrix(0, nrow = nrow(bin_mat), ncol = ncol(bin_mat))
  
  # MATLAB: prob_s = (1:length(p_simu))'/length(p_simu);
  # length(p_simu) in MATLAB here is number of rows (simulations), because rows >> cols.
  n_simu <- nrow(p_simu)
  prob_s <- seq_len(n_simu) / n_simu
  
  for (i in seq_len(width_s)) {
    in_val <- sort(p_simu[, i], decreasing = FALSE)
    cum_p_bins <- cumsum(bin_mat[, i])
    
    # MATLAB interp1(prob_s, in_val, cum_p_bins)
    p_bins[, i] <- approx(x = prob_s, y = in_val, xout = cum_p_bins, rule = 2)$y
  }
  
  p_bins
}
