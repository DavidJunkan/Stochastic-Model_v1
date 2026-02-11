#install.packages("moments")
library(readxl)
library(moments)

n_simu    <- 15000
y_base    <- 2024
mu        <- 0
val_miss  <- -9.999
n_spline  <- 1000

root_path <- "D:/Dropbox/Stochastic Models/Code_replication/R_version"
data_path <- "D:/Dropbox/Stochastic Models/Code_replication/R_version/Data"
setwd(root_path)
source("gen_logn.R")
source("f_detrend.R")
source("gen_pbin.R")
source("ukernel_gauss.R")

## Reading the Data, including Deterministic Case and Historical Data 
det_input <- "deterministic_input.xlsx"
T_read    <- read_excel(file.path(data_path, det_input),sheet = 1)
crop_name <- colnames(T_read)[-1]
T_read[is.na(T_read)] <- val_miss
p_det_dat <- T_read
y_p_det   <- as.matrix(T_read[[1]])
p_det     <- as.matrix(T_read[, -1])
n_crop    <- ncol(p_det)
n_yp_det  <- length(y_p_det)

T_read    <- read_excel(file.path(data_path, det_input),sheet = 2)
v_det_dat <- T_read
sigma     <- as.matrix(T_read[1, -1])
sigma[is.na(sigma)]       <- 0
time_adj  <- as.matrix(T_read[2, -1])
time_adj[is.na(time_adj)] <- 1

T_read    <- read_excel(file.path(data_path, det_input), sheet = 3)
T_read[is.na(T_read)] <- val_miss
a_det_dat <- T_read
y_a_det   <- as.matrix(T_read[[1]])                 
acre_det  <- as.matrix(T_read[, -1])   

T_read    <- read_excel(file.path(data_path, det_input), sheet = 4)
q_det_dat <- T_read
T_read[is.na(T_read)] <- val_miss
y_q_det   <- as.matrix(T_read[[1]])
prod_det  <- as.matrix(T_read[, -1])

T_read    <- read_excel(file.path(data_path, "bins.xlsx"), sheet = 1)
bins      <- as.matrix(T_read[[2]])
bins      <- bins[1:(length(bins) - 2)]

seeds     <- read.csv(file.path(data_path, "seedmat3.csv"))

## Simulate the Price Deviation
p_val     <- p_det[1, ]
seed_val  <- seeds[1, 1]
res       <- gen_logn(mu, sigma, time_adj, p_val, n_simu, seed_val)
p_simu_dv <- res$p_dev
p_simu    <- res$p_simu
check_p   <- colMeans(p_simu) - p_val

## Reading All Production and Acre Data
n_year    <- 50
year_q    <- (y_base - n_year + 1):y_base

acre_dat  <- matrix(0, nrow = n_year, ncol = n_crop)
prod_dat  <- matrix(0, nrow = n_year, ncol = n_crop)

for (i in seq_len(n_crop)) {
  
  # Read sheet i from acre_input.xlsx and keep columns 2:4
  T_acre <- read_excel(file.path(data_path, "acre_input.xlsx"), sheet = i)
  A_acre <- T_acre[, 2:4]
  A_acre[] <- lapply(A_acre, as.numeric)   # ensure numeric comparisons work
  
  # Read sheet i from production_input.xlsx and keep columns 2:4
  T_prod <- read_excel(file.path(data_path, "production_input.xlsx"), sheet = i)
  A_prod <- T_prod[, 2:4]
  A_prod[] <- lapply(A_prod, as.numeric)
  
  for (j in seq_len(n_year)) {
    
    # MATLAB: find(A(:,1)==year & A(:,2)==99)
    set_1 <- which(A_acre[[1]] == year_q[j] & A_acre[[2]] == 99)
    set_2 <- which(A_prod[[1]] == year_q[j] & A_prod[[2]] == 99)
    
    if (length(set_1) > 0) {
      acre_dat[j, i] <- A_acre[[3]][set_1[1]]
    } else {
      acre_dat[j, i] <- val_miss
    }
    
    if (length(set_2) > 0) {
      prod_dat[j, i] <- A_prod[[3]][set_2[1]]
    } else {
      prod_dat[j, i] <- val_miss
    }
  }
}


##  Read in price deviation data
T_read      <- read_excel(file.path(data_path, det_input), sheet = 5)
T_read[is.na(T_read)] <- val_miss
year_p      <- T_read[[1]]
p_dv_dat    <- as.matrix(T_read[, -1])
colnames(p_dv_dat) <- paste0(crop_name, "_price")


##  Yield Simulation, Kernel
yld_val     <- prod_dat / acre_dat
dt_method   <- 2
res         <- f_detrend(yld_val, dt_method)

yld_dt      <- res$yld_dt
yld_t       <- res$yld_t
yld_dv      <- res$yld_dv
colnames(yld_dv) <- paste0(crop_name, "_yield")

yld_simu    <- matrix(0, nrow = n_simu, ncol = n_crop)
yld_simu_dv <- matrix(0, nrow = n_simu, ncol = n_crop)

for (i in seq_len(n_crop)) {
  
  # Use detrended series for crop i
  z <- yld_dt[, i]
  # Recommended: drop placeholder missing values before KDE
  z <- z[is.finite(z) & z != val_miss]
  yld_min <- min(z)
  yld_max <- max(z)
  ext_gap <- 2 * sd(z)
  yld_ssp <- seq(from = yld_min - ext_gap, to = yld_max + ext_gap, length.out = n_spline)
  # Kernel density on the grid using your Gaussian kernel function
  kres    <- ukernel_gauss(x = yld_ssp, z = z, h = 0, w = rep(1, length(z)))
  pdf     <- kres$f
  # Build CDF by rectangle rule
  dx      <- (max(yld_ssp) - min(yld_ssp)) / (n_spline - 1)
  P       <- cumsum(pdf * dx)
  # Draw uniforms and invert CDF via linear interpolation
  U       <- runif(n_simu)
  y       <- approx(x = P, y = yld_ssp, xout = U, rule = 2)$y
  # Clamp to bounds
  y       <- pmax(y, min(yld_ssp))
  y       <- pmin(y, max(yld_ssp))
  yld_simu[, i]    <- y
  m       <- mean(yld_simu[, i])
  yld_simu_dv[, i] <- (yld_simu[, i] - m) / m
}

# Column-wise moments: REAL data
mean_real  <- colMeans(yld_dt)
std_real   <- apply(yld_dt, 2, sd)
kurt_real  <- apply(yld_dt, 2, kurtosis)
skew_real  <- apply(yld_dt, 2, skewness)

# Column-wise moments: SIMULATED data
mean_simu  <- colMeans(yld_simu)
std_simu   <- apply(yld_simu, 2, sd)
kurt_simu  <- apply(yld_simu, 2, kurtosis)
skew_simu  <- apply(yld_simu, 2, skewness)

check_simu <- rbind(mean_real,std_real,kurt_real,skew_real,
                    mean_simu,std_simu,kurt_simu,skew_simu)
colnames(check_simu) <- crop_name
print(check_simu)


## Capula Method
year_cor   <- sort(intersect(year_p, year_q))
id_p       <- match(year_cor, year_p)   # row indices in year_p
id_q       <- match(year_cor, year_q)   # row indices in year_q
p_dv_cor   <- p_dv_dat[id_p, , drop = FALSE]
yld_dv_cor <- yld_dv[id_q, , drop = FALSE]

data_cor   <- cbind(p_dv_cor[, 1:(n_crop - 2), drop = FALSE], yld_dv_cor)
rho        <- cor(data_cor)
pp         <- chol(rho)  # upper-triangular in R, same convention as MATLAB
n_cor      <- ncol(rho)
XX_cor     <- matrix(rnorm(n_simu * n_cor), nrow = n_simu, ncol = n_cor)
RXX_cor    <- XX_cor %*% pp
UXX_cor    <- pnorm(RXX_cor, mean = 0, sd = 1)

simu_raw   <- cbind(p_simu_dv[, 1:(n_crop - 2), drop = FALSE],yld_simu_dv)
simu_cor   <- matrix(0, nrow = nrow(simu_raw), ncol = ncol(simu_raw))
prob_ref   <- seq_len(n_simu) / n_simu

for (i in seq_len(n_cor)) {
  
  x_val <- sort(simu_raw[, i], decreasing = FALSE)
  
  # MATLAB uses prob_ref = (1:n_simu)/n_simu, and clamps small U to first interval.
  # rule = 2 handles boundary behavior; we also enforce the same lower bound as MATLAB.
  u     <- UXX_cor[, i]
  u[u <= 1 / n_simu] <- 1 / n_simu
  simu_cor[, i] <- approx(x = prob_ref, y = x_val, xout = u, rule = 2)$y
  
  # center the column (same as simu_cor(:,i) = simu_cor(:,i) - mean(...))
  simu_cor[, i] <- simu_cor[, i] - mean(simu_cor[, i])
}

rho_simu  <- cor(simu_cor)
check_rho <- rho_simu - rho
print(check_rho)

p_simu_dv_cor   <- simu_cor[, 1:(n_crop - 2), drop = FALSE]
yld_simu_dv_cor <- simu_cor[, (n_crop - 1):ncol(data_cor), drop = FALSE]


## Now use the baseline deterministic values
#  You are free to decide which year as your target
y_target   <- 2024
id_p       <- which(y_p_det == y_target)
id_a       <- which(y_a_det == y_target)
id_q       <- which(y_q_det == y_target)

p_target   <- p_det[id_p[1], , drop = FALSE]
arc_target <- acre_det[id_a[1], , drop = FALSE]
yld_target <- prod_det[id_q[1], , drop = FALSE] / arc_target

# --- Allocate ---
yld_det_simu_iid  <- matrix(0, nrow = n_simu, ncol = n_crop)
yld_det_simu_cor  <- matrix(0, nrow = n_simu, ncol = n_crop)
prod_det_simu_iid <- matrix(0, nrow = n_simu, ncol = n_crop)
prod_det_simu_cor <- matrix(0, nrow = n_simu, ncol = n_crop)
p_det_simu_iid    <- matrix(0, nrow = n_simu, ncol = n_crop)
p_det_simu_cor    <- matrix(0, nrow = n_simu, ncol = n_crop)

# --- Main loop ---
for (i in seq_len(n_crop)) {
  
  yld_det_simu_iid[, i]  <- yld_simu_dv[, i]      * yld_target[1, i] + yld_target[1, i]
  yld_det_simu_cor[, i]  <- yld_simu_dv_cor[, i]  * yld_target[1, i] + yld_target[1, i]
  prod_det_simu_iid[, i] <- yld_det_simu_iid[, i] * arc_target[1, i]
  prod_det_simu_cor[, i] <- yld_det_simu_cor[, i] * arc_target[1, i]
  
  #  i_crop = 6 is rice
  id_p_dv <- if (i <= (n_crop - 2)) i else 6
  p_det_simu_iid[, i]    <- p_simu_dv[, i] * p_target[1, i] + p_target[1, i]
  p_det_simu_cor[, i]    <- p_simu_dv_cor[, id_p_dv] * p_target[1, i] + p_target[1, i]
}

# --- Column-wise moments (MATLAB-equivalent) ---
mean_simu_iid <- c(colMeans(p_det_simu_iid), colMeans(yld_det_simu_iid))
std_simu_iid  <- c(apply(p_det_simu_iid, 2, sd), apply(yld_det_simu_iid, 2, sd))
kurt_simu_iid <- c(apply(p_det_simu_iid, 2, kurtosis), apply(yld_det_simu_iid, 2, kurtosis))
skew_simu_iid <- c(apply(p_det_simu_iid, 2, skewness), apply(yld_det_simu_iid, 2, skewness))

mean_simu_cor <- c(colMeans(p_det_simu_cor), colMeans(yld_det_simu_cor))
std_simu_cor  <- c(apply(p_det_simu_cor, 2, sd), apply(yld_det_simu_cor, 2, sd))
kurt_simu_cor <- c(apply(p_det_simu_cor, 2, kurtosis), apply(yld_det_simu_cor, 2, kurtosis))
skew_simu_cor <- c(apply(p_det_simu_cor, 2, skewness), apply(yld_det_simu_cor, 2, skewness))

check_cor <- rbind(mean_simu_iid, std_simu_iid, kurt_simu_iid, skew_simu_iid,
                   mean_simu_cor,std_simu_cor,kurt_simu_cor,skew_simu_cor)

print(check_cor)


## Generate distribution by bins (same target year)
idx_bin <- 1

bins_p_cor   <- gen_pbin(p_det_simu_cor,  bins, idx_bin)
bins_yld_cor <- gen_pbin(yld_det_simu_cor, bins, idx_bin)
