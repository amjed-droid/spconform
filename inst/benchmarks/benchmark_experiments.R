# Benchmark Experiment: Big Data Acceleration in spconform
# Comparing:
# 1. Original (recomputing bandwidth in loop)
# 2. Optimized Baseline (Precomputed bandwidth + Vectorized)
# 3. Proposal 1: k-NN Local Calibration (k = 50, 100, 250)
# 4. Proposal 2: Truncated / Compact Support Kernel
# 5. Proposal 3: Multi-Core Parallel Processing

source("R/scp_geostatistical.R")
source("R/utils.R")
source("R/methods.R")

# 1. Optimized Baseline (Precomputed bandwidth once)
scp_geostatistical_opt_base <- function(s_train, y_train, s0, pred_fun,
                                        alpha = 0.1, split = 0.5,
                                        bandwidth = NULL, seed = NULL) {
  s_train <- as.matrix(s_train)
  s0 <- as.matrix(s0)
  n <- nrow(s_train)
  m <- nrow(s0)
  
  if (!is.null(seed)) set.seed(seed)
  n_fit <- floor(split * n)
  fit_idx <- sample(seq_len(n), n_fit)
  cal_idx <- setdiff(seq_len(n), fit_idx)
  
  s_fit <- s_train[fit_idx, , drop = FALSE]
  y_fit <- y_train[fit_idx]
  s_cal <- s_train[cal_idx, , drop = FALSE]
  y_cal <- y_train[cal_idx]
  n_cal <- length(cal_idx)
  
  # Precompute bandwidth ONCE (using smart subsampling if n_cal is large)
  if (is.null(bandwidth)) {
    sub_sample <- if (n_cal > 1000) s_cal[sample(n_cal, 500), , drop = FALSE] else s_cal
    bandwidth <- stats::median(stats::dist(sub_sample))
    if (!is.finite(bandwidth) || bandwidth <= 0) bandwidth <- 1
  }
  h2_2 <- 2 * bandwidth^2
  
  cal_pred <- pred_fun(s_fit, y_fit, s_cal)
  resid_abs <- abs(y_cal - cal_pred)
  pred0 <- pred_fun(s_fit, y_fit, s0)
  
  ord <- order(resid_abs)
  r_sorted <- resid_abs[ord]
  
  lower <- upper <- numeric(m)
  q_level <- min(1, (1 - alpha) * (n_cal + 1) / n_cal)
  
  for (j in seq_len(m)) {
    # Distance to calibration set
    d2 <- (s_cal[, 1] - s0[j, 1])^2 + (s_cal[, 2] - s0[j, 2])^2
    w <- exp(-d2 / h2_2) + 1e-12
    w_sorted <- w[ord]
    cw <- cumsum(w_sorted) / sum(w_sorted)
    
    q_idx <- which(cw >= q_level)[1]
    if (is.na(q_idx)) q_idx <- n_cal
    qhat <- r_sorted[q_idx]
    
    lower[j] <- pred0[j] - qhat
    upper[j] <- pred0[j] + qhat
  }
  
  list(pred = pred0, lower = lower, upper = upper, alpha = alpha)
}

# 2. Proposal 1: Fast k-NN Local Calibration
scp_geostatistical_knn <- function(s_train, y_train, s0, pred_fun,
                                   alpha = 0.1, split = 0.5,
                                   k_neighbors = 100, bandwidth = NULL, seed = NULL) {
  s_train <- as.matrix(s_train)
  s0 <- as.matrix(s0)
  n <- nrow(s_train)
  m <- nrow(s0)
  
  if (!is.null(seed)) set.seed(seed)
  n_fit <- floor(split * n)
  fit_idx <- sample(seq_len(n), n_fit)
  cal_idx <- setdiff(seq_len(n), fit_idx)
  
  s_fit <- s_train[fit_idx, , drop = FALSE]
  y_fit <- y_train[fit_idx]
  s_cal <- s_train[cal_idx, , drop = FALSE]
  y_cal <- y_train[cal_idx]
  n_cal <- length(cal_idx)
  k <- min(k_neighbors, n_cal)
  
  cal_pred <- pred_fun(s_fit, y_fit, s_cal)
  resid_abs <- abs(y_cal - cal_pred)
  pred0 <- pred_fun(s_fit, y_fit, s0)
  
  # Local bandwidth
  if (is.null(bandwidth)) {
    sub_sample <- if (n_cal > 1000) s_cal[sample(n_cal, 500), , drop = FALSE] else s_cal
    bandwidth <- stats::median(stats::dist(sub_sample))
    if (!is.finite(bandwidth) || bandwidth <= 0) bandwidth <- 1
  }
  h2_2 <- 2 * bandwidth^2
  
  lower <- upper <- numeric(m)
  q_level <- min(1, (1 - alpha) * (k + 1) / k)
  
  for (j in seq_len(m)) {
    d2 <- (s_cal[, 1] - s0[j, 1])^2 + (s_cal[, 2] - s0[j, 2])^2
    knn_idx <- order(d2)[seq_len(k)]
    
    sub_d2 <- d2[knn_idx]
    sub_r <- resid_abs[knn_idx]
    
    w <- exp(-sub_d2 / h2_2) + 1e-12
    ord <- order(sub_r)
    r_sorted <- sub_r[ord]
    w_sorted <- w[ord]
    cw <- cumsum(w_sorted) / sum(w_sorted)
    
    q_idx <- which(cw >= q_level)[1]
    if (is.na(q_idx)) q_idx <- k
    qhat <- r_sorted[q_idx]
    
    lower[j] <- pred0[j] - qhat
    upper[j] <- pred0[j] + qhat
  }
  
  list(pred = pred0, lower = lower, upper = upper, alpha = alpha)
}

# 3. Proposal 2: Truncated Support Kernel (Cutoff at 3*bandwidth)
scp_geostatistical_truncated <- function(s_train, y_train, s0, pred_fun,
                                         alpha = 0.1, split = 0.5,
                                         bandwidth = NULL, seed = NULL) {
  s_train <- as.matrix(s_train)
  s0 <- as.matrix(s0)
  n <- nrow(s_train)
  m <- nrow(s0)
  
  if (!is.null(seed)) set.seed(seed)
  n_fit <- floor(split * n)
  fit_idx <- sample(seq_len(n), n_fit)
  cal_idx <- setdiff(seq_len(n), fit_idx)
  
  s_fit <- s_train[fit_idx, , drop = FALSE]
  y_fit <- y_train[fit_idx]
  s_cal <- s_train[cal_idx, , drop = FALSE]
  y_cal <- y_train[cal_idx]
  n_cal <- length(cal_idx)
  
  cal_pred <- pred_fun(s_fit, y_fit, s_cal)
  resid_abs <- abs(y_cal - cal_pred)
  pred0 <- pred_fun(s_fit, y_fit, s0)
  
  if (is.null(bandwidth)) {
    sub_sample <- if (n_cal > 1000) s_cal[sample(n_cal, 500), , drop = FALSE] else s_cal
    bandwidth <- stats::median(stats::dist(sub_sample))
    if (!is.finite(bandwidth) || bandwidth <= 0) bandwidth <- 1
  }
  cutoff_d2 <- (3 * bandwidth)^2
  h2_2 <- 2 * bandwidth^2
  
  lower <- upper <- numeric(m)
  
  for (j in seq_len(m)) {
    d2 <- (s_cal[, 1] - s0[j, 1])^2 + (s_cal[, 2] - s0[j, 2])^2
    in_range <- which(d2 <= cutoff_d2)
    
    if (length(in_range) < 10) {
      in_range <- order(d2)[1:min(20, n_cal)]
    }
    
    sub_d2 <- d2[in_range]
    sub_r <- resid_abs[in_range]
    k_local <- length(in_range)
    
    w <- exp(-sub_d2 / h2_2) + 1e-12
    ord <- order(sub_r)
    r_sorted <- sub_r[ord]
    w_sorted <- w[ord]
    cw <- cumsum(w_sorted) / sum(w_sorted)
    
    q_level <- min(1, (1 - alpha) * (k_local + 1) / k_local)
    q_idx <- which(cw >= q_level)[1]
    if (is.na(q_idx)) q_idx <- k_local
    qhat <- r_sorted[q_idx]
    
    lower[j] <- pred0[j] - qhat
    upper[j] <- pred0[j] + qhat
  }
  
  list(pred = pred0, lower = lower, upper = upper, alpha = alpha)
}

cat("=================================================================\n")
cat("          SPCONFORM BIG DATA BENCHMARKING EXPERIMENT             \n")
cat("=================================================================\n\n")

# Dataset: 5,000 points train, 2,000 test locations
set.seed(42)
N_train <- 5000
M_test <- 2000

cat(sprintf("Training Sample N: %d | Test Locations M: %d\n\n", N_train, M_test))

s_train <- matrix(runif(2 * N_train, min = 0, max = 10), ncol = 2)
true_surface <- function(s) sin(s[, 1] / 1.5) + cos(s[, 2] / 1.5) + 0.3 * (s[, 1] * s[, 2]) / 20
y_train <- true_surface(s_train) + rnorm(N_train, sd = 0.4)

s0 <- matrix(runif(2 * M_test, min = 0, max = 10), ncol = 2)
y0_true <- true_surface(s0) + rnorm(M_test, sd = 0.4)

pred_fun <- function(s_tr, y_tr, s_nw) {
  fit <- lm(y_tr ~ s_tr[, 1] + s_tr[, 2] + I(s_tr[, 1]^2) + I(s_tr[, 2]^2))
  as.numeric(cbind(1, s_nw[, 1], s_nw[, 2], s_nw[, 1]^2, s_nw[, 2]^2) %*% coef(fit))
}

# 1. Test Optimized Baseline
cat("1. Testing Optimized Full-Sample Baseline...\n")
t0 <- proc.time()
res_base <- scp_geostatistical_opt_base(s_train, y_train, s0, pred_fun, alpha = 0.1, seed = 123)
t_base <- (proc.time() - t0)["elapsed"]
cov_base <- mean(y0_true >= res_base$lower & y0_true <= res_base$upper)
width_base <- mean(res_base$upper - res_base$lower)

# 2. Test k-NN (k = 50)
cat("2. Testing Proposal 1: k-NN (k = 50)...\n")
t0 <- proc.time()
res_knn50 <- scp_geostatistical_knn(s_train, y_train, s0, pred_fun, alpha = 0.1, k_neighbors = 50, seed = 123)
t_knn50 <- (proc.time() - t0)["elapsed"]
cov_knn50 <- mean(y0_true >= res_knn50$lower & y0_true <= res_knn50$upper)
width_knn50 <- mean(res_knn50$upper - res_knn50$lower)

# 3. Test k-NN (k = 100)
cat("3. Testing Proposal 1: k-NN (k = 100)...\n")
t0 <- proc.time()
res_knn100 <- scp_geostatistical_knn(s_train, y_train, s0, pred_fun, alpha = 0.1, k_neighbors = 100, seed = 123)
t_knn100 <- (proc.time() - t0)["elapsed"]
cov_knn100 <- mean(y0_true >= res_knn100$lower & y0_true <= res_knn100$upper)
width_knn100 <- mean(res_knn100$upper - res_knn100$lower)

# 4. Test k-NN (k = 250)
cat("4. Testing Proposal 1: k-NN (k = 250)...\n")
t0 <- proc.time()
res_knn250 <- scp_geostatistical_knn(s_train, y_train, s0, pred_fun, alpha = 0.1, k_neighbors = 250, seed = 123)
t_knn250 <- (proc.time() - t0)["elapsed"]
cov_knn250 <- mean(y0_true >= res_knn250$lower & y0_true <= res_knn250$upper)
width_knn250 <- mean(res_knn250$upper - res_knn250$lower)

# 5. Test Truncated Kernel
cat("5. Testing Proposal 2: Compact / Truncated Support Kernel...\n")
t0 <- proc.time()
res_trunc <- scp_geostatistical_truncated(s_train, y_train, s0, pred_fun, alpha = 0.1, seed = 123)
t_trunc <- (proc.time() - t0)["elapsed"]
cov_trunc <- mean(y0_true >= res_trunc$lower & y0_true <= res_trunc$upper)
width_trunc <- mean(res_trunc$upper - res_trunc$lower)

cat("\n=================================================================\n")
cat("                       EXPERIMENTAL RESULTS                      \n")
cat("=================================================================\n")

res_df <- data.frame(
  Method = c("Optimized Full Calibration (2500 pts)", 
             "k-NN (k = 50)", 
             "k-NN (k = 100)", 
             "k-NN (k = 250)", 
             "Truncated Kernel (3 * Bandwidth)"),
  Time_Sec = round(c(t_base, t_knn50, t_knn100, t_knn250, t_trunc), 3),
  Speedup = sprintf("%.1fx", c(1.0, t_base / t_knn50, t_base / t_knn100, t_base / t_knn250, t_base / t_trunc)),
  Empirical_Coverage = sprintf("%.2f%%", c(cov_base, cov_knn50, cov_knn100, cov_knn250, cov_trunc) * 100),
  Mean_Interval_Width = round(c(width_base, width_knn50, width_knn100, width_knn250, width_trunc), 4)
)
print(res_df, row.names = FALSE)
cat("=================================================================\n")
