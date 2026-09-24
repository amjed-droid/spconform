pkgload::load_all(".")

cat("=================================================================\n")
cat("   MEGA-STRESS TEST: 20,000 TRAINING & 5,000 TARGET LOCATIONS     \n")
cat("=================================================================\n\n")

set.seed(42)
N_train <- 20000
M_target <- 5000

cat(sprintf("Simulating massive dataset: N = %d, M = %d...\n", N_train, M_target))

s_train <- matrix(runif(2 * N_train, min = 0, max = 100), ncol = 2)
colnames(s_train) <- c("x", "y")
f_spatial <- function(s) sin(s[, 1] / 10) * cos(s[, 2] / 10)
y_train <- f_spatial(s_train) + rnorm(N_train, sd = 0.25)

s0_test <- matrix(runif(2 * M_target, min = 0, max = 100), ncol = 2)
colnames(s0_test) <- c("x", "y")
y_test_true <- f_spatial(s0_test) + rnorm(M_target, sd = 0.25)

pred_fun <- function(s_tr, y_tr, s_new) {
  fit <- lm(y_tr ~ s_tr[, 1] + s_tr[, 2])
  as.numeric(cbind(1, s_new[, 1], s_new[, 2]) %*% coef(fit))
}

cat("Running localized spatial calibration (k = 50, n_cores = 4)...\n")
t0 <- Sys.time()
res_mega <- scp_geostatistical(
  s_train = s_train,
  y_train = y_train,
  s0 = s0_test,
  pred_fun = pred_fun,
  alpha = 0.10,
  k_neighbors = 50,
  n_cores = 4L,
  progress = TRUE,
  seed = 123
)
t_elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
rep_mega <- coverage_report(res_mega, y_test_true)

cat(sprintf("\nTotal Execution Time: %.3f seconds\n", t_elapsed))
cat(sprintf("Achieved Coverage:    %.2f%% (Nominal: 90.00%%)\n", rep_mega$coverage * 100))
cat(sprintf("Mean Interval Width:  %.4f\n", rep_mega$mean_width))
cat("Sample Tidy Predictions Output:\n")
print(head(as.data.frame(res_mega), 3))
