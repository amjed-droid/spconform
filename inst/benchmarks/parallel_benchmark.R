# Test Parallel Execution across CPU Cores
source("R/scp_geostatistical.R")
source("R/utils.R")
source("R/methods.R")
source("tests/benchmark_experiments.R")

cat("\n=================================================================\n")
cat("       PARALLEL MULTI-CORE TEST (Windows & Unix Compatible)      \n")
cat("=================================================================\n\n")

num_cores <- parallel::detectCores(logical = FALSE)
cat(sprintf("Physical CPU Cores Detected: %d cores\n\n", num_cores))

# Parallel version of k-NN
scp_geostatistical_parallel <- function(s_train, y_train, s0, pred_fun,
                                        alpha = 0.1, split = 0.5,
                                        k_neighbors = 100, n_cores = 2, seed = 123) {
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
  
  sub_sample <- if (n_cal > 1000) s_cal[sample(n_cal, 500), , drop = FALSE] else s_cal
  bandwidth <- stats::median(stats::dist(sub_sample))
  if (!is.finite(bandwidth) || bandwidth <= 0) bandwidth <- 1
  h2_2 <- 2 * bandwidth^2
  
  q_level <- min(1, (1 - alpha) * (k + 1) / k)
  
  # Setup cluster
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl))
  
  # Export needed data
  parallel::clusterExport(cl, varlist = c("s_cal", "s0", "resid_abs", "h2_2", "k", "q_level"), envir = environment())
  
  # Worker function
  calc_point <- function(j) {
    d2 <- (s_cal[, 1] - s0[j, 1])^2 + (s_cal[, 2] - s0[j, 2])^2
    knn_idx <- order(d2)[seq_len(k)]
    sub_d2 <- d2[knn_idx]
    sub_r <- resid_abs[knn_idx]
    w <- exp(-sub_d2 / h2_2) + 1e-12
    ord <- order(sub_r)
    cw <- cumsum(w[ord]) / sum(w)
    q_idx <- which(cw >= q_level)[1]
    if (is.na(q_idx)) q_idx <- k
    sub_r[ord][q_idx]
  }
  
  qhats <- unlist(parallel::parLapply(cl, seq_len(m), calc_point))
  
  list(pred = pred0, lower = pred0 - qhats, upper = pred0 + qhats, alpha = alpha)
}

# Test on 50,000 points
set.seed(42)
N_train <- 50000
M_test <- 5000

s_train <- matrix(runif(2 * N_train, min = 0, max = 20), ncol = 2)
true_surface <- function(s) sin(s[, 1] / 2) + cos(s[, 2] / 2)
y_train <- true_surface(s_train) + rnorm(N_train, sd = 0.5)

s0 <- matrix(runif(2 * M_test, min = 0, max = 20), ncol = 2)
y0_true <- true_surface(s0) + rnorm(M_test, sd = 0.5)

pred_fun <- function(s_tr, y_tr, s_nw) {
  fit <- lm(y_tr ~ s_tr[, 1] + s_tr[, 2])
  as.numeric(cbind(1, s_nw[, 1], s_nw[, 2]) %*% coef(fit))
}

cat("Running Single-Core (1 Core)...\n")
t0 <- proc.time()
res_1core <- scp_geostatistical_knn(s_train, y_train, s0, pred_fun, alpha = 0.1, k_neighbors = 100, seed = 123)
t_1core <- (proc.time() - t0)["elapsed"]
cat(sprintf(">> Single-Core Time: %.2f sec\n", t_1core))

cat(sprintf("Running Parallel on (%d Cores)...\n", num_cores))
t0 <- proc.time()
res_multi <- scp_geostatistical_parallel(s_train, y_train, s0, pred_fun, alpha = 0.1, k_neighbors = 100, n_cores = num_cores, seed = 123)
t_multi <- (proc.time() - t0)["elapsed"]
cat(sprintf(">> Multi-Core Time: %.2f sec (Speedup: %.1fx)\n\n", t_multi, t_1core / t_multi))

cat("=================================================================\n")
