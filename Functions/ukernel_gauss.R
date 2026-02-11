ukernel_gauss <- function(x, z, h, w) {
  
  x <- as.numeric(x)
  z <- as.numeric(z)
  w <- as.numeric(w)
  
  n_x <- length(x)
  n_z <- length(z)
  
  f <- numeric(n_x)
  d <- numeric(n_x)
  
  # Calculate the bandwidth (only if h <= 0)
  if (h <= 0) {
    s <- sqrt(stats::var(z))
    
    q1 <- as.numeric(stats::quantile(z, probs = 0.25, type = 7, na.rm = TRUE))
    q3 <- as.numeric(stats::quantile(z, probs = 0.75, type = 7, na.rm = TRUE))
    iqr <- q3 - q1
    
    h <- 0.9 * min(s, iqr / 1.34) / (n_z^0.2)
  }
  
  for (i in seq_len(n_x)) {
    arg <- (x[i] - z) / h
    
    kff <- stats::dnorm(arg)          # N(0,1) pdf
    kfd <- arg * stats::dnorm(arg)    # arg .* normpdf(arg)
    
    f[i] <- mean(kff * w) / h
    d[i] <- mean(kfd * w) / (h^2)
  }
  
  list(f = f, d = d, h = h)
}
