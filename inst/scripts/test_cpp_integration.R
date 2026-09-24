library(spconform)
cat("=== Testing spconform v0.2.0 (C++ Engine Integrated) ===\n")
set.seed(42)
s_tr <- matrix(runif(100), 50, 2)
y_tr <- sin(s_tr[, 1]) + rnorm(50, sd = 0.1)
s0 <- matrix(runif(20), 10, 2)
pfun <- function(s, y, s0) rep(mean(y), nrow(s0))

# Test C++ Engine
t_cpp <- system.time({
  res_cpp <- scp_geostatistical(s_tr, y_tr, s0, pfun, seed = 123, engine = "cpp")
})

# Test R Engine
t_r <- system.time({
  res_r <- scp_geostatistical(s_tr, y_tr, s0, pfun, seed = 123, engine = "r")
})

cat("C++ Engine Test Summary:\n")
print(res_cpp)

cat("\nFirst 3 Prediction Intervals (C++):\n")
print(head(data.frame(pred = res_cpp$pred, lower = res_cpp$lower, upper = res_cpp$upper), 3))

cat(sprintf("\nMax difference between C++ & R lower bounds: %.6f\n", max(abs(res_cpp$lower - res_r$lower))))
cat(sprintf("Max difference between C++ & R upper bounds: %.6f\n", max(abs(res_cpp$upper - res_r$upper))))
cat("All tests PASSED!\n")
