f_detrend <- function(yld, method) {
  
  yld <- as.matrix(yld)
  n_crop <- ncol(yld)
  n_t <- nrow(yld)
  
  yld_t  <- matrix(0, nrow = n_t, ncol = n_crop)  # trend
  yld_c  <- matrix(0, nrow = n_t, ncol = n_crop)  # cycle/residual
  yld_dv <- matrix(0, nrow = n_t, ncol = n_crop)  # deviation ratio
  yld_dt <- matrix(0, nrow = n_t, ncol = n_crop)  # detrended (scaled)
  
  for (i in seq_len(n_crop)) {
    
    if (method == 1) {
      v_max <- max(yld[, i], na.rm = TRUE)
      v_min <- min(yld[, i], na.rm = TRUE)
      
      yld_t[, i] <- seq(from = v_min, to = v_max, length.out = n_t)
      yld_c[, i] <- yld[, i] - yld_t[, i]
      
    } else if (method == 2) {
      T <- seq_len(n_t)
      X <- cbind(1, T)  # [ones, T]
      
      # MATLAB: fitlm(X, y, 'Intercept', false) -> OLS with no extra intercept
      reg <- lm(yld[, i] ~ X - 1)
      
      yld_t[, i] <- as.numeric(X %*% coef(reg))
      yld_c[, i] <- yld[, i] - yld_t[, i]
      
    } else if (method == 3) {
      # HP filter: trend + cycle
      # Choose lambda appropriate to your data frequency (see note below)
      if (!requireNamespace("mFilter", quietly = TRUE)) {
        stop("Package 'mFilter' is required for method == 3. Install via install.packages('mFilter').")
      }
      hp <- mFilter::hpfilter(yld[, i], freq = 100)  # freq=lambda; adjust if needed
      
      yld_t[, i] <- as.numeric(hp$trend)
      yld_c[, i] <- as.numeric(hp$cycle)
      
    } else {
      stop("method must be 1, 2, or 3.")
    }
    
    yld_dv[, i] <- yld_c[, i] / yld_t[, i]
    yld_dt[, i] <- (yld[, i] / yld_t[, i]) * yld_t[n_t, i]
  }
  
  list(yld_dt = yld_dt, yld_t = yld_t, yld_dv = yld_dv)
}
