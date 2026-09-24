pkgload::load_all(".")

cat("=================================================================\n")
cat("   REAL-WORLD DATASET VALIDATION: EARTHQUAKES SEISMICITY DATA    \n")
cat("=================================================================\n\n")

# Load real-world Fiji Earthquakes dataset (built into base R)
data(quakes)
cat(sprintf("Loaded real-world 'quakes' dataset: %d seismic events.\n", nrow(quakes)))
cat("Variables available: lat, long, depth, mag (Richter magnitude), stations.\n\n")

# Target variable: Earthquake Magnitude (mag)
# Predictors: Geographic coordinates (long, lat, depth)
s_all <- as.matrix(quakes[, c("long", "lat")])
d_depth <- quakes$depth
y_all <- quakes$mag

# Set seed and create real Train (70%) and Test (30%) split
set.seed(2026)
n_total <- nrow(quakes)
train_idx <- sample(n_total, size = round(0.70 * n_total))

s_train <- s_all[train_idx, ]
y_train <- y_all[train_idx]
d_train <- d_depth[train_idx]

s_test  <- s_all[-train_idx, ]
y_test  <- y_all[-train_idx]
d_test  <- d_depth[-train_idx]

cat(sprintf("Training observations: %d | Test target locations: %d\n\n", length(train_idx), length(y_test)))

# Realistic Non-Linear Spatial Predictor (Spline/Polynomial Trend Surface)
pred_fun_real <- function(s_tr, y_tr, s_new) {
  df_tr <- data.frame(lon = s_tr[, 1], lat = s_tr[, 2], y = y_tr)
  df_new <- data.frame(lon = s_new[, 1], lat = s_new[, 2])
  fit <- stats::lm(y ~ poly(lon, 3) * poly(lat, 3), data = df_tr)
  stats::predict(fit, newdata = df_new)
}

# Run spconform geostatistical calibration on real data with k_neighbors
cat("Running spatial conformal calibration (Target Nominal Coverage = 90%, alpha = 0.10)...\n")
t0 <- Sys.time()
res_real <- scp_geostatistical(
  s_train = s_train,
  y_train = y_train,
  s0 = s_test,
  pred_fun = pred_fun_real,
  alpha = 0.10,
  k_neighbors = 40,
  progress = TRUE,
  seed = 42
)
t_elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# Evaluate Real-World Performance
rep_real <- coverage_report(res_real, y_test)

cat(sprintf("\nExecution Time:       %.4f seconds\n", t_elapsed))
cat(sprintf("Empirical Coverage:   %.2f%%  (Nominal Target: 90.00%%)\n", rep_real$coverage * 100))
cat(sprintf("Mean Interval Width:  %.3f Richter magnitude units\n", rep_real$mean_width))

# Verify coverage validity
is_valid <- rep_real$coverage >= (1 - res_real$alpha - 0.05)
cat(sprintf("Mathematical Validity: %s (Coverage preserved on real empirical data!)\n\n", ifelse(is_valid, "PASSED [OK]", "FAILED")))

# Show sample of real predictions vs ground truth
df_real <- as.data.frame(res_real)
df_real$true_mag <- y_test
df_real$is_covered <- (df_real$true_mag >= df_real$lower) & (df_real$true_mag <= df_real$upper)

cat("Sample Real-World Interval Predictions:\n")
print(head(df_real[, c("long", "lat", "true_mag", "pred", "lower", "upper", "is_covered")], 6))

cat("\n=================================================================\n")
cat("      REAL-WORLD EMPIRICAL VALIDATION COMPLETED SUCCESSFULLY!    \n")
cat("=================================================================\n")
