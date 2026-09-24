# Massive Dataset Scale Test: N = 50,000 | M = 5,000
source("R/scp_geostatistical.R")
source("R/utils.R")
source("R/methods.R")

source("tests/benchmark_experiments.R") # loads helper functions

cat("\n=================================================================\n")
cat("          MASSIVE SCALE TEST: N = 50,000 | M = 5,000             \n")
cat("=================================================================\n\n")

set.seed(42)
N_train <- 50000
M_test <- 5000

cat("Generating 50,000 spatial observations...\n")
s_train <- matrix(runif(2 * N_train, min = 0, max = 20), ncol = 2)
true_surface <- function(s) sin(s[, 1] / 2) + cos(s[, 2] / 2)
y_train <- true_surface(s_train) + rnorm(N_train, sd = 0.5)

s0 <- matrix(runif(2 * M_test, min = 0, max = 20), ncol = 2)
y0_true <- true_surface(s0) + rnorm(M_test, sd = 0.5)

pred_fun <- function(s_tr, y_tr, s_nw) {
  fit <- lm(y_tr ~ s_tr[, 1] + s_tr[, 2])
  as.numeric(cbind(1, s_nw[, 1], s_nw[, 2]) %*% coef(fit))
}

# 1. Test k-NN (k = 100)
cat("Testing Fast k-NN (k = 100) on 50,000 data points...\n")
t0 <- proc.time()
res_knn <- scp_geostatistical_knn(s_train, y_train, s0, pred_fun, alpha = 0.1, k_neighbors = 100, seed = 123)
t_knn <- (proc.time() - t0)["elapsed"]
cov_knn <- mean(y0_true >= res_knn$lower & y0_true <= res_knn$upper)
width_knn <- mean(res_knn$upper - res_knn$lower)

cat(sprintf(">> k-NN (k=100) Finished in: %.2f sec | Coverage: %.2f%% | Mean Width: %.3f\n", 
            t_knn, cov_knn * 100, width_knn))

# 2. Test Optimized Full Calibration (25,000 calibration points per target)
cat("\nTesting Full Calibration (All 25,000 points per target)...\n")
t0 <- proc.time()
res_base <- scp_geostatistical_opt_base(s_train, y_train, s0, pred_fun, alpha = 0.1, seed = 123)
t_base <- (proc.time() - t0)["elapsed"]
cov_base <- mean(y0_true >= res_base$lower & y0_true <= res_base$upper)
width_base <- mean(res_base$upper - res_base$lower)

cat(sprintf(">> Full Calibration Finished in: %.2f sec | Coverage: %.2f%% | Mean Width: %.3f\n\n", 
            t_base, cov_base * 100, width_base))

cat("=================================================================\n")
cat(sprintf("SPEEDUP ON 50,000 POINTS: %.1fx FASTER!\n", t_base / t_knn))
cat("=================================================================\n")
