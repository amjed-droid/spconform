## =============================================================================
## Script: eval_areal_heldout.R
## Description: Held-out areal evaluation and Monte Carlo cross-validation
##              for localized conformal prediction via scp_areal().
##
## Author: Ahmed Sattar Jabbar
## Package: spconform (v0.1.0)
## License: GPL (>= 3)
## =============================================================================

library(spconform)
library(sp)

## -----------------------------------------------------------------------------
## 1. DATA PREPARATION & AREAL GRID AGGREGATION
## -----------------------------------------------------------------------------

## Load Meuse river benchmark data
data(meuse, package = "sp")

## Aggregate continuous coordinates onto a regular 6x6 spatial lattice
xbreaks <- seq(min(meuse$x), max(meuse$x), length.out = 7)
ybreaks <- seq(min(meuse$y), max(meuse$y), length.out = 7)

meuse$cell_x  <- cut(meuse$x, xbreaks, include.lowest = TRUE, labels = FALSE)
meuse$cell_y  <- cut(meuse$y, ybreaks, include.lowest = TRUE, labels = FALSE)
meuse$cell_id <- (meuse$cell_y - 1) * 6 + meuse$cell_x

## Compute mean log-zinc concentration for each occupied grid cell
agg <- stats::aggregate(log(zinc) ~ cell_id, data = meuse, FUN = mean)
names(agg) <- c("cell_id", "y")

## Merge spatial cell coordinates and sort by cell ID
cell_coords <- unique(meuse[, c("cell_id", "cell_x", "cell_y")])
agg <- merge(agg, cell_coords, by = "cell_id")
agg <- agg[order(agg$cell_id), ]

n_cells <- nrow(agg)

## Build full binary adjacency matrix (Queen's contiguity rule)
adj_full <- matrix(0, nrow = n_cells, ncol = n_cells)
for (i in 1:n_cells) {
  for (j in 1:n_cells) {
    if (i != j) {
      dx <- abs(agg$cell_x[i] - agg$cell_x[j])
      dy <- abs(agg$cell_y[i] - agg$cell_y[j])
      if (dx <= 1 && dy <= 1) {
        adj_full[i, j] <- 1
      }
    }
  }
}

## -----------------------------------------------------------------------------
## 2. MONTE CARLO EVALUATION (50 RANDOM TRAIN / TEST SPLITS)
## -----------------------------------------------------------------------------

n_reps <- 50
alpha  <- 0.2  # Nominal 80% coverage target

results <- data.frame(
  rep            = seq_len(n_reps),
  train_coverage = numeric(n_reps),
  train_width    = numeric(n_reps),
  test_coverage  = numeric(n_reps),
  test_width     = numeric(n_reps),
  n_train        = integer(n_reps),
  n_test         = integer(n_reps)
)

for (r in seq_len(n_reps)) {
  set.seed(r)
  
  ## 70/30 random partition of areal units
  train_idx <- sample(n_cells, size = floor(0.7 * n_cells))
  test_idx  <- setdiff(seq_len(n_cells), train_idx)
  
  y_train <- agg$y[train_idx]
  y_test  <- agg$y[test_idx]
  
  adj_train <- adj_full[train_idx, train_idx]
  
  ## (A) Training set LOO calibration via scp_areal
  cal_out <- tryCatch(
    scp_areal(y_train, adjacency = adj_train, alpha = alpha, decay = 0.5),
    error = function(e) NULL
  )
  
  if (is.null(cal_out)) next
  
  ## Training LOO empirical coverage & mean interval width
  train_covered <- (cal_out$lower <= y_train) & (y_train <= cal_out$upper)
  train_covered <- train_covered[!is.na(train_covered)]
  
  if (length(train_covered) > 0) {
    results$train_coverage[r] <- mean(train_covered)
    results$train_width[r]    <- mean(cal_out$upper - cal_out$lower, na.rm = TRUE)
  }
  
  ## Extract leave-one-out nonconformity absolute residuals
  m <- length(train_idx)
  cal_scores <- numeric(m)
  for (i in seq_along(train_idx)) {
    idx_loo <- setdiff(seq_along(train_idx), i)
    if (length(idx_loo) > 0) {
      pred_loo <- mean(y_train[idx_loo])
      cal_scores[i] <- abs(y_train[i] - pred_loo)
    } else {
      cal_scores[i] <- 0
    }
  }
  
  ## (B) Held-out test unit prediction & localized conformal calibration
  tau <- min(1, (1 - alpha) * (m + 1) / m)
  
  test_lower <- numeric(length(test_idx))
  test_upper <- numeric(length(test_idx))
  test_pred  <- numeric(length(test_idx))
  
  for (j in seq_along(test_idx)) {
    test_i <- test_idx[j]
    
    ## Breadth-First Search (BFS) for shortest graph hops on full adjacency
    visited <- rep(FALSE, n_cells)
    dist    <- rep(Inf, n_cells)
    queue   <- test_i
    dist[test_i]    <- 0
    visited[test_i] <- TRUE
    
    qhead <- 1
    while (qhead <= length(queue)) {
      current <- queue[qhead]
      qhead   <- qhead + 1
      neighbors <- which(adj_full[current, ] == 1)
      for (nb in neighbors) {
        if (!visited[nb]) {
          visited[nb] <- TRUE
          dist[nb]    <- dist[current] + 1
          queue       <- c(queue, nb)
        }
      }
    }
    
    ## Graph-distance exponential decay kernel weights
    train_dist <- dist[train_idx]
    weights    <- exp(-0.5 * train_dist)
    
    ## Point prediction: local neighbor mean (fallback to global training mean)
    adj_to_train <- adj_full[test_i, train_idx]
    if (sum(adj_to_train) > 0) {
      test_pred[j] <- mean(y_train[adj_to_train == 1])
    } else {
      test_pred[j] <- mean(y_train)
    }
    
    ## Weighted conformal empirical quantile calculation
    ord            <- order(cal_scores)
    sorted_scores  <- cal_scores[ord]
    sorted_weights <- weights[ord]
    
    if (sum(sorted_weights) > 0) {
      cum_weights <- cumsum(sorted_weights) / sum(sorted_weights)
      q_idx       <- min(which(cum_weights >= tau))
      q_hat       <- sorted_scores[q_idx]
    } else {
      q_hat       <- max(cal_scores, na.rm = TRUE)
    }
    
    test_lower[j] <- test_pred[j] - q_hat
    test_upper[j] <- test_pred[j] + q_hat
  }
  
  ## Compute empirical coverage on held-out test units
  test_covered             <- (y_test >= test_lower) & (y_test <= test_upper)
  results$test_coverage[r] <- mean(test_covered, na.rm = TRUE)
  results$test_width[r]    <- mean(test_upper - test_lower, na.rm = TRUE)
  results$n_train[r]       <- length(train_idx)
  results$n_test[r]        <- length(test_idx)
}

## -----------------------------------------------------------------------------
## 3. SUMMARY OUTPUT & DIAGNOSTICS
## -----------------------------------------------------------------------------

cat("\n=========================================================\n")
cat("   Held-out Areal Evaluation Summary (50 Replications)   \n")
cat("=========================================================\n\n")

cat("Training Set (LOO Calibration Baseline):\n")
cat(sprintf("  - Mean Coverage: %.3f (SD: %.3f)\n", 
            mean(results$train_coverage, na.rm = TRUE), 
            sd(results$train_coverage, na.rm = TRUE)))
cat(sprintf("  - Mean Width:    %.3f (SD: %.3f)\n\n", 
            mean(results$train_width, na.rm = TRUE), 
            sd(results$train_width, na.rm = TRUE)))

cat("Held-out Test Units (Out-of-sample Generalization):\n")
cat(sprintf("  - Mean Coverage: %.3f (SD: %.3f)\n", 
            mean(results$test_coverage, na.rm = TRUE), 
            sd(results$test_coverage, na.rm = TRUE)))
cat(sprintf("  - Mean Width:    %.3f (SD: %.3f)\n\n", 
            mean(results$test_width, na.rm = TRUE), 
            sd(results$test_width, na.rm = TRUE)))

cat(sprintf("Completed Splits: %d / %d\n", 
            sum(!is.na(results$test_coverage)), n_reps))
cat("=========================================================\n\n")

## -----------------------------------------------------------------------------
## 4. VISUALIZATION & DIAGNOSTIC PLOTS
## -----------------------------------------------------------------------------

oldpar <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1.5))

## Panel A: Empirical Coverage Distribution
boxplot(
  list("Training (LOO)" = results$train_coverage, 
       "Held-out Test"  = results$test_coverage),
  main = "Coverage: Training LOO vs. Held-out Test",
  ylab = "Empirical Coverage",
  col  = c("#A6CEE3", "#B2DF8A"),
  ylim = c(0.4, 1.0),
  las  = 1
)
abline(h = 1 - alpha, col = "red", lty = 2, lwd = 2)
legend("bottomright", legend = sprintf("Nominal (%.2f)", 1 - alpha), 
       col = "red", lty = 2, lwd = 2, bty = "n")

## Panel B: Prediction Interval Width Distribution
boxplot(
  list("Training (LOO)" = results$train_width,
       "Held-out Test"  = results$test_width),
  main = "Interval Width: Training vs. Test",
  ylab = "Mean Interval Width",
  col  = c("#A6CEE3", "#B2DF8A"),
  las  = 1
)

par(oldpar)