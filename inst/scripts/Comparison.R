## =============================================================
## spconform - Model Comparison Script
## Compares coverage and interval width for different base predictors
## =============================================================

library(sp)
library(spconform)
library(mgcv)      # For GAM
library(ranger)    # For Random Forest

## -------------------------------------------------------------
## 1. Load real data: Meuse river dataset
## -------------------------------------------------------------
data(meuse)

s <- as.matrix(meuse[, c("x", "y")])
y <- log(meuse$zinc)

## -------------------------------------------------------------
## 2. Define base predictors
## -------------------------------------------------------------

## 2.1. Linear model (quadratic trend) - Misspecified / Weak
pred_fun_lm <- function(s_train, y_train, s_new) {
  fit <- lm(y_train ~ s_train[, 1] + s_train[, 2] +
              I(s_train[, 1]^2) + I(s_train[, 2]^2))
  cbind(1, s_new[, 1], s_new[, 2],
        s_new[, 1]^2, s_new[, 2]^2) %*% coef(fit)
}

## 2.2. GAM (Generalized Additive Model) - Better fit
pred_fun_gam <- function(s_train, y_train, s_new) {
  train_df <- data.frame(x = s_train[, 1], y = s_train[, 2], z = y_train)
  fit <- gam(z ~ s(x) + s(y), data = train_df)  # Smooth terms for coordinates
  new_df <- data.frame(x = s_new[, 1], y = s_new[, 2])
  as.numeric(predict(fit, newdata = new_df))
}

## 2.3. Random Forest - Flexible / Strong
pred_fun_rf <- function(s_train, y_train, s_new) {
  train_df <- data.frame(x = s_train[, 1], y = s_train[, 2], z = y_train)
  fit <- ranger(z ~ x + y, data = train_df, num.trees = 500, 
                mtry = 1, min.node.size = 5, seed = 42)
  new_df <- data.frame(x = s_new[, 1], y = s_new[, 2])
  predict(fit, data = new_df)$predictions
}

## -------------------------------------------------------------
## 3. Perform comparison across random splits
## -------------------------------------------------------------
set.seed(123)
n_splits <- 100  # Increase to 50 for more stability
n <- nrow(s)

# Storage for results
results_lm <- data.frame(coverage = numeric(n_splits), width = numeric(n_splits))
results_gam <- data.frame(coverage = numeric(n_splits), width = numeric(n_splits))
results_rf  <- data.frame(coverage = numeric(n_splits), width = numeric(n_splits))

for (i in 1:n_splits) {
  # Split data (70/30)
  idx <- sample(n, floor(0.7 * n))
  s_train <- s[idx, ]; y_train <- y[idx]
  s_test  <- s[-idx, ]; y_test  <- y[-idx]
  
  # LM
  out_lm <- scp_geostatistical(s_train, y_train, s_test, pred_fun_lm, 
                               alpha = 0.1, seed = i)
  rep_lm <- coverage_report(out_lm, y_test)
  results_lm$coverage[i] <- rep_lm$coverage
  results_lm$width[i] <- rep_lm$mean_width
  
  # GAM
  out_gam <- scp_geostatistical(s_train, y_train, s_test, pred_fun_gam, 
                                alpha = 0.1, seed = i)
  rep_gam <- coverage_report(out_gam, y_test)
  results_gam$coverage[i] <- rep_gam$coverage
  results_gam$width[i] <- rep_gam$mean_width
  
  # Random Forest
  out_rf <- scp_geostatistical(s_train, y_train, s_test, pred_fun_rf, 
                               alpha = 0.1, seed = i)
  rep_rf <- coverage_report(out_rf, y_test)
  results_rf$coverage[i] <- rep_rf$coverage
  results_rf$width[i] <- rep_rf$mean_width
}

## -------------------------------------------------------------
## 4. Print summary statistics
## -------------------------------------------------------------
cat("\n===== Comparison of Base Predictors (Target coverage: 90%) =====\n")

cat("\n--- Linear Model (Quadratic Trend) ---\n")
cat(sprintf("Mean coverage: %.3f (SD: %.3f)\n", 
            mean(results_lm$coverage), sd(results_lm$coverage)))
cat(sprintf("Mean interval width: %.3f\n", mean(results_lm$width)))

cat("\n--- GAM (Smooth terms) ---\n")
cat(sprintf("Mean coverage: %.3f (SD: %.3f)\n", 
            mean(results_gam$coverage), sd(results_gam$coverage)))
cat(sprintf("Mean interval width: %.3f\n", mean(results_gam$width)))

cat("\n--- Random Forest ---\n")
cat(sprintf("Mean coverage: %.3f (SD: %.3f)\n", 
            mean(results_rf$coverage), sd(results_rf$coverage)))
cat(sprintf("Mean interval width: %.3f\n", mean(results_rf$width)))

## -------------------------------------------------------------
## 5. Visualize the results
## -------------------------------------------------------------

# Boxplot comparison for coverage
boxplot(list(LM = results_lm$coverage, 
             GAM = results_gam$coverage, 
             RF = results_rf$coverage),
        main = "Empirical Coverage Comparison (Target = 0.90)",
        ylab = "Coverage", col = c("lightblue", "lightgreen", "lightpink"))
abline(h = 0.90, col = "red", lty = 2, lwd = 2)
legend("bottomright", legend = "Nominal target (0.90)", col = "red", lty = 2, bty = "n")

# Boxplot comparison for width
boxplot(list(LM = results_lm$width, 
             GAM = results_gam$width, 
             RF = results_rf$width),
        main = "Interval Width Comparison",
        ylab = "Mean Interval Width", 
        col = c("lightblue", "lightgreen", "lightpink"))