

``` r
## =============================================================================
## spconform - Master Reproducible Analysis Script 
## Generates all figures (fig1.pdf - fig8.pdf), tables, and benchmark outputs 
## =============================================================================

options(stringsAsFactors = FALSE)

SEED <- 123
set.seed(SEED)

OUTPUT_DIR <- "figures"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

library(sp)
library(spconform)
library(mgcv)       # For Spatio-temporal GAMs
library(ranger)     # For Random Forest benchmarks
library(bmstdr)     # For New York Ozone dataset
```

``` r
## PART 1: MEUSE RIVER DATA - GEOSTATISTICAL ILLUSTRATION
```

``` r
cat("\n============================================================\n")
```

```
## 
## ============================================================
```

``` r
cat("  PART 1: Meuse River Data - Geostatistical Analysis        \n")
```

```
##   PART 1: Meuse River Data - Geostatistical Analysis
```

``` r
cat("============================================================\n")
```

```
## ============================================================
```

``` r
## 1.1 Load data and define base predictor
data(meuse, package = "sp")
s <- as.matrix(meuse[, c("x", "y")])
y <- log(meuse$zinc)

pred_fun <- function(s_train, y_train, s_new) {
  fit <- lm(y_train ~ s_train[, 1] + s_train[, 2] +
              I(s_train[, 1]^2) + I(s_train[, 2]^2))
  cbind(1, s_new[, 1], s_new[, 2],
        s_new[, 1]^2, s_new[, 2]^2) %*% coef(fit)
}

## 1.2 Figure 1: Spatial layout
pdf(file.path(OUTPUT_DIR, "fig1.pdf"), width = 6, height = 5)
plot(meuse$x, meuse$y,
     col = rgb(0.2, 0.4, 0.8, 0.5), pch = 19,
     xlab = "X coordinate", ylab = "Y coordinate",
     main = "Meuse River Sampling Locations")
dev.off()
```

```
## RStudioGD 
##         2
```

``` r
## 1.3 Single 70/30 split
set.seed(SEED)
n <- nrow(s)
idx <- sample(n, floor(0.7 * n))
s_train <- s[idx, ]; y_train <- y[idx]
s_test  <- s[-idx, ]; y_test  <- y[-idx]

out <- scp_geostatistical(s_train, y_train, s_test, pred_fun,
                          alpha = 0.1, seed = SEED)

## Figure 2: Single-split prediction intervals
pdf(file.path(OUTPUT_DIR, "fig2.pdf"), width = 7, height = 5)
if (any(is.na(out$lower)) || any(is.na(out$upper))) {
  valid <- !is.na(out$lower) & !is.na(out$upper)
  out_plot <- out
  out_plot$lower <- out$lower[valid]
  out_plot$upper <- out$upper[valid]
  out_plot$pred  <- out$pred[valid]
  plot(out_plot, y_true = y_test[valid])
} else {
  plot(out, y_true = y_test)
}
dev.off()
```

```
## RStudioGD 
##         2
```

``` r
## 1.4 Monte Carlo: 50 random splits
set.seed(SEED)
coverages <- numeric(50)
widths    <- numeric(50)

for (i in 1:50) {
  idx_i  <- sample(n, floor(0.7 * n))
  s_tr   <- s[idx_i, ]; y_tr <- y[idx_i]
  s_te   <- s[-idx_i, ]; y_te <- y[-idx_i]
  out_i  <- scp_geostatistical(s_tr, y_tr, s_te, pred_fun,
                               alpha = 0.1, seed = i)
  rep_i  <- coverage_report(out_i, y_te)
  coverages[i] <- rep_i$coverage
  widths[i]    <- rep_i$mean_width
}

## Figure 3: Coverage histogram across 50 splits
pdf(file.path(OUTPUT_DIR, "fig3.pdf"), width = 6, height = 5)
hist(coverages, breaks = 15, col = "lightblue", border = "white",
     main = "Empirical Coverage Across 50 Random Splits",
     xlab = "Empirical Coverage", xlim = c(0.7, 1))
abline(v = 0.90, col = "red", lwd = 2, lty = 2)
legend("topleft", legend = "Nominal target (0.90)",
       col = "red", lty = 2, bty = "n")
dev.off()
```

```
## RStudioGD 
##         2
```

``` r
## 1.5 Figure 4: Spatial distribution of interval width
plot_df <- data.frame(x = s_test[, 1], y = s_test[, 2],
                      width = out$upper - out$lower)
pdf(file.path(OUTPUT_DIR, "fig4.pdf"), width = 6, height = 5)
plot(plot_df$x, plot_df$y, cex = plot_df$width, pch = 19,
     col = rgb(0.2, 0.4, 0.8, 0.5),
     xlab = "X coordinate", ylab = "Y coordinate",
     main = "Spatial Distribution of Interval Width")
dev.off()
```

```
## RStudioGD 
##         2
```

``` r
## 1.6 Summary statistics for Geostatistical section
cat(sprintf("Single split coverage: %.3f\n", coverage_report(out, y_test)$coverage))
```

```
## Single split coverage: 1.000
```

``` r
cat(sprintf("Mean coverage (50 splits): %.3f (SD: %.3f)\n", 
            mean(coverages), sd(coverages)))
```

```
## Mean coverage (50 splits): 0.920 (SD: 0.044)
```

``` r
cat(sprintf("Mean width (50 splits): %.3f\n", mean(widths)))
```

```
## Mean width (50 splits): 1.973
```

``` r
## PART 2: DIAGNOSTIC REPORT (FIGURE 5)
```

``` r
cat("\n============================================================\n")
```

```
## 
## ============================================================
```

``` r
cat("  PART 2: Spatial Diagnostics (Figure 5)                    \n")
```

```
##   PART 2: Spatial Diagnostics (Figure 5)
```

``` r
cat("============================================================\n")
```

```
## ============================================================
```

``` r
pdf(file.path(OUTPUT_DIR, "fig8.pdf"), width = 8.5, height = 7)
diag_meuse <- diagnose(
  object  = out,
  y_true  = y_test,
  s_test  = s_test,
  n_bins  = 4,
  plot    = TRUE
)
```

```
## Note: Empirical coverage (1) exceeds nominal (0.9) by >5%. Consider reducing 'bandwidth' for tighter intervals.
```

``` r
dev.off()
```

```
## RStudioGD 
##         2
```

``` r
saveRDS(diag_meuse, file = "spconform_diagnostics.rds")
cat("Figure 5 (diagnostics) saved to:", file.path(OUTPUT_DIR, "fig8.pdf"), "\n")
```

```
## Figure 5 (diagnostics) saved to: figures/fig8.pdf
```

``` r
## PART 3: MEUSE RIVER DATA - AREAL ILLUSTRATION
```

``` r
cat("\n============================================================\n")
```

```
## 
## ============================================================
```

``` r
cat("  PART 3: Meuse River Data - Areal Lattice Analysis         \n")
```

```
##   PART 3: Meuse River Data - Areal Lattice Analysis
```

``` r
cat("============================================================\n")
```

```
## ============================================================
```

``` r
## 3.1 Aggregate Meuse to 6x6 grid
xbreaks <- seq(min(meuse$x), max(meuse$x), length.out = 7)
ybreaks <- seq(min(meuse$y), max(meuse$y), length.out = 7)

meuse$cell_x  <- cut(meuse$x, xbreaks, include.lowest = TRUE, labels = FALSE)
meuse$cell_y  <- cut(meuse$y, ybreaks, include.lowest = TRUE, labels = FALSE)
meuse$cell_id <- (meuse$cell_y - 1) * 6 + meuse$cell_x

agg <- aggregate(log(zinc) ~ cell_id, data = meuse, FUN = mean)
names(agg) <- c("cell_id", "y")

cell_coords <- unique(meuse[, c("cell_id", "cell_x", "cell_y")])
agg <- merge(agg, cell_coords, by = "cell_id")
agg <- agg[order(agg$cell_id), ]

n_cells <- nrow(agg)
adj <- matrix(0, n_cells, n_cells)
for (i in 1:n_cells) {
  for (j in 1:n_cells) {
    if (i != j) {
      dx <- abs(agg$cell_x[i] - agg$cell_x[j])
      dy <- abs(agg$cell_y[i] - agg$cell_y[j])
      if (dx <= 1 && dy <= 1) adj[i, j] <- 1
    }
  }
}

## 3.2 Areal conformal prediction
out2 <- scp_areal(agg$y, adjacency = adj, alpha = 0.2, decay = 0.5)

## Figure 6 (in paper): Areal prediction intervals
pdf(file.path(OUTPUT_DIR, "fig5.pdf"), width = 7, height = 5)
if (any(is.na(out2$lower)) || any(is.na(out2$upper))) {
  valid <- !is.na(out2$lower) & !is.na(out2$upper)
  out2_plot <- out2
  out2_plot$lower <- out2$lower[valid]
  out2_plot$upper <- out2$upper[valid]
  out2_plot$pred  <- out2$pred[valid]
  plot(out2_plot, y_true = agg$y[valid])
} else {
  plot(out2, y_true = agg$y)
}
dev.off()
```

```
## RStudioGD 
##         2
```

``` r
## Figure 7 (in paper): Interval width comparison (Geostatistical vs Areal)
pdf(file.path(OUTPUT_DIR, "fig6.pdf"), width = 6, height = 5)
geo_width   <- out$upper - out$lower
areal_width <- out2$upper - out2$lower
geo_width   <- geo_width[!is.na(geo_width)]
areal_width <- areal_width[!is.na(areal_width)]

if (length(geo_width) > 0 && length(areal_width) > 0) {
  boxplot(list(Geostatistical = geo_width,
               Areal           = areal_width),
          main = "Interval Width Comparison",
          ylab = "Interval Width",
          col  = c("lightblue", "lightgreen"))
}
dev.off()
```

```
## RStudioGD 
##         2
```

``` r
## Summary statistics for Areal section
if (any(is.na(out2$lower)) || any(is.na(out2$upper))) {
  valid <- !is.na(out2$lower) & !is.na(out2$upper)
  rep_areal <- coverage_report(
    list(pred = out2$pred[valid], lower = out2$lower[valid], 
         upper = out2$upper[valid], alpha = out2$alpha),
    agg$y[valid]
  )
} else {
  rep_areal <- coverage_report(out2, agg$y)
}
cat(sprintf("Areal coverage: %.3f\n", rep_areal$coverage))
```

```
## Areal coverage: 0.810
```

``` r
cat(sprintf("Areal mean width: %.3f\n", rep_areal$mean_width))
```

```
## Areal mean width: 1.790
```

``` r
## PART 4: HELD-OUT AREAL EVALUATION (FIGURE 8)
```

``` r
cat("\n============================================================\n")
```

```
## 
## ============================================================
```

``` r
cat("  PART 4: Held-out Areal Evaluation (Figure 8)              \n")
```

```
##   PART 4: Held-out Areal Evaluation (Figure 8)
```

``` r
cat("============================================================\n")
```

```
## ============================================================
```

``` r
script_path <- file.path("inst", "scripts", "eval_areal_heldout.R")
if (!file.exists(script_path)) {
  script_path <- system.file("scripts", "eval_areal_heldout.R", package = "spconform")
}

if (file.exists(script_path)) {
  pdf(file.path(OUTPUT_DIR, "fig7.pdf"), width = 8.5, height = 4.5)
  source(script_path, local = TRUE)
  dev.off()
  cat("Figure 8 saved to:", file.path(OUTPUT_DIR, "fig7.pdf"), "\n")
} else {
  warning("eval_areal_heldout.R not found in expected paths!")
}
```

```
## Warning in scp_areal(y_train, adjacency = adj_train, alpha = alpha, decay = 0.5): Unit 8 has negligible
## neighbourhood weights; falling back to unweighted quantile.
```

```
## Warning in scp_areal(y_train, adjacency = adj_train, alpha = alpha, decay = 0.5): Unit 4 has negligible
## neighbourhood weights; falling back to unweighted quantile.
```

```
## Warning in scp_areal(y_train, adjacency = adj_train, alpha = alpha, decay = 0.5): Unit 10 has negligible
## neighbourhood weights; falling back to unweighted quantile.
```

```
## Warning in scp_areal(y_train, adjacency = adj_train, alpha = alpha, decay = 0.5): Unit 7 has negligible
## neighbourhood weights; falling back to unweighted quantile.
```

```
## Warning in scp_areal(y_train, adjacency = adj_train, alpha = alpha, decay = 0.5): Unit 6 has negligible
## neighbourhood weights; falling back to unweighted quantile.
```

```
## Warning in scp_areal(y_train, adjacency = adj_train, alpha = alpha, decay = 0.5): Unit 13 has negligible
## neighbourhood weights; falling back to unweighted quantile.
## Warning in scp_areal(y_train, adjacency = adj_train, alpha = alpha, decay = 0.5): Unit 13 has negligible
## neighbourhood weights; falling back to unweighted quantile.
```

```
## Warning in scp_areal(y_train, adjacency = adj_train, alpha = alpha, decay = 0.5): Unit 3 has negligible
## neighbourhood weights; falling back to unweighted quantile.
```

```
## 
## =========================================================
##    Held-out Areal Evaluation Summary (50 Replications)   
## =========================================================
## 
## Training Set (LOO Calibration Baseline):
##   - Mean Coverage: 0.834 (SD: 0.042)
##   - Mean Width:    2.075 (SD: 0.467)
## 
## Held-out Test Units (Out-of-sample Generalization):
##   - Mean Coverage: 0.857 (SD: 0.135)
##   - Mean Width:    1.974 (SD: 0.258)
## 
## Completed Splits: 50 / 50
## =========================================================
```

```
## Figure 8 saved to: figures/fig7.pdf
```

``` r
## PART 5: SENSITIVITY ANALYSIS: MODEL COMPARISON (100 SPLITS)
```

``` r
cat("\n============================================================\n")
```

```
## 
## ============================================================
```

``` r
cat("  PART 5: Sensitivity Analysis (Linear vs GAM vs RF)        \n")
```

```
##   PART 5: Sensitivity Analysis (Linear vs GAM vs RF)
```

``` r
cat("============================================================\n")
```

```
## ============================================================
```

``` r
pred_fun_lm <- function(s_train, y_train, s_new) {
  fit <- lm(y_train ~ s_train[, 1] + s_train[, 2] +
              I(s_train[, 1]^2) + I(s_train[, 2]^2))
  cbind(1, s_new[, 1], s_new[, 2],
        s_new[, 1]^2, s_new[, 2]^2) %*% coef(fit)
}

pred_fun_gam <- function(s_train, y_train, s_new) {
  train_df <- data.frame(x = s_train[, 1], y = s_train[, 2], z = y_train)
  fit <- gam(z ~ s(x) + s(y), data = train_df)
  new_df <- data.frame(x = s_new[, 1], y = s_new[, 2])
  as.numeric(predict(fit, newdata = new_df))
}

pred_fun_rf <- function(s_train, y_train, s_new) {
  train_df <- data.frame(x = s_train[, 1], y = s_train[, 2], z = y_train)
  fit <- ranger(z ~ x + y, data = train_df, num.trees = 500,
                mtry = 1, min.node.size = 5, seed = SEED)
  new_df <- data.frame(x = s_new[, 1], y = s_new[, 2])
  predict(fit, data = new_df)$predictions
}

set.seed(SEED)
n_splits <- 100

results_lm  <- data.frame(coverage = numeric(n_splits), width = numeric(n_splits))
results_gam <- data.frame(coverage = numeric(n_splits), width = numeric(n_splits))
results_rf  <- data.frame(coverage = numeric(n_splits), width = numeric(n_splits))

pb <- txtProgressBar(min = 0, max = n_splits, style = 3)
```

```
##   |                                                                                                                   |                                                                                                           |   0%
```

``` r
for (i in 1:n_splits) {
  idx <- sample(n, floor(0.7 * n))
  s_train <- s[idx, ]; y_train <- y[idx]
  s_test  <- s[-idx, ]; y_test  <- y[-idx]
  
  out_lm <- scp_geostatistical(s_train, y_train, s_test, pred_fun_lm, alpha = 0.1, seed = i)
  rep_lm <- coverage_report(out_lm, y_test)
  results_lm$coverage[i] <- rep_lm$coverage
  results_lm$width[i]    <- rep_lm$mean_width
  
  out_gam <- scp_geostatistical(s_train, y_train, s_test, pred_fun_gam, alpha = 0.1, seed = i)
  rep_gam <- coverage_report(out_gam, y_test)
  results_gam$coverage[i] <- rep_gam$coverage
  results_gam$width[i]    <- rep_gam$mean_width
  
  out_rf <- scp_geostatistical(s_train, y_train, s_test, pred_fun_rf, alpha = 0.1, seed = i)
  rep_rf <- coverage_report(out_rf, y_test)
  results_rf$coverage[i] <- rep_rf$coverage
  results_rf$width[i]    <- rep_rf$mean_width
  
  setTxtProgressBar(pb, i)
}
```

```
##   |                                                                                                                   |=                                                                                                          |   1%  |                                                                                                                   |==                                                                                                         |   2%  |                                                                                                                   |===                                                                                                        |   3%  |                                                                                                                   |====                                                                                                       |   4%  |                                                                                                                   |=====                                                                                                      |   5%  |                                                                                                                   |======                                                                                                     |   6%  |                                                                                                                   |=======                                                                                                    |   7%  |                                                                                                                   |=========                                                                                                  |   8%  |                                                                                                                   |==========                                                                                                 |   9%  |                                                                                                                   |===========                                                                                                |  10%  |                                                                                                                   |============                                                                                               |  11%  |                                                                                                                   |=============                                                                                              |  12%  |                                                                                                                   |==============                                                                                             |  13%  |                                                                                                                   |===============                                                                                            |  14%  |                                                                                                                   |================                                                                                           |  15%  |                                                                                                                   |=================                                                                                          |  16%  |                                                                                                                   |==================                                                                                         |  17%  |                                                                                                                   |===================                                                                                        |  18%  |                                                                                                                   |====================                                                                                       |  19%  |                                                                                                                   |=====================                                                                                      |  20%  |                                                                                                                   |======================                                                                                     |  21%  |                                                                                                                   |========================                                                                                   |  22%  |                                                                                                                   |=========================                                                                                  |  23%  |                                                                                                                   |==========================                                                                                 |  24%  |                                                                                                                   |===========================                                                                                |  25%  |                                                                                                                   |============================                                                                               |  26%  |                                                                                                                   |=============================                                                                              |  27%  |                                                                                                                   |==============================                                                                             |  28%  |                                                                                                                   |===============================                                                                            |  29%  |                                                                                                                   |================================                                                                           |  30%  |                                                                                                                   |=================================                                                                          |  31%  |                                                                                                                   |==================================                                                                         |  32%  |                                                                                                                   |===================================                                                                        |  33%  |                                                                                                                   |====================================                                                                       |  34%  |                                                                                                                   |=====================================                                                                      |  35%  |                                                                                                                   |=======================================                                                                    |  36%  |                                                                                                                   |========================================                                                                   |  37%  |                                                                                                                   |=========================================                                                                  |  38%  |                                                                                                                   |==========================================                                                                 |  39%  |                                                                                                                   |===========================================                                                                |  40%  |                                                                                                                   |============================================                                                               |  41%  |                                                                                                                   |=============================================                                                              |  42%  |                                                                                                                   |==============================================                                                             |  43%  |                                                                                                                   |===============================================                                                            |  44%  |                                                                                                                   |================================================                                                           |  45%  |                                                                                                                   |=================================================                                                          |  46%  |                                                                                                                   |==================================================                                                         |  47%  |                                                                                                                   |===================================================                                                        |  48%  |                                                                                                                   |====================================================                                                       |  49%  |                                                                                                                   |======================================================                                                     |  50%  |                                                                                                                   |=======================================================                                                    |  51%  |                                                                                                                   |========================================================                                                   |  52%  |                                                                                                                   |=========================================================                                                  |  53%  |                                                                                                                   |==========================================================                                                 |  54%  |                                                                                                                   |===========================================================                                                |  55%  |                                                                                                                   |============================================================                                               |  56%  |                                                                                                                   |=============================================================                                              |  57%  |                                                                                                                   |==============================================================                                             |  58%  |                                                                                                                   |===============================================================                                            |  59%  |                                                                                                                   |================================================================                                           |  60%  |                                                                                                                   |=================================================================                                          |  61%  |                                                                                                                   |==================================================================                                         |  62%  |                                                                                                                   |===================================================================                                        |  63%  |                                                                                                                   |====================================================================                                       |  64%  |                                                                                                                   |======================================================================                                     |  65%  |                                                                                                                   |=======================================================================                                    |  66%  |                                                                                                                   |========================================================================                                   |  67%  |                                                                                                                   |=========================================================================                                  |  68%  |                                                                                                                   |==========================================================================                                 |  69%  |                                                                                                                   |===========================================================================                                |  70%  |                                                                                                                   |============================================================================                               |  71%  |                                                                                                                   |=============================================================================                              |  72%  |                                                                                                                   |==============================================================================                             |  73%  |                                                                                                                   |===============================================================================                            |  74%  |                                                                                                                   |================================================================================                           |  75%  |                                                                                                                   |=================================================================================                          |  76%  |                                                                                                                   |==================================================================================                         |  77%  |                                                                                                                   |===================================================================================                        |  78%  |                                                                                                                   |=====================================================================================                      |  79%  |                                                                                                                   |======================================================================================                     |  80%  |                                                                                                                   |=======================================================================================                    |  81%  |                                                                                                                   |========================================================================================                   |  82%  |                                                                                                                   |=========================================================================================                  |  83%  |                                                                                                                   |==========================================================================================                 |  84%  |                                                                                                                   |===========================================================================================                |  85%  |                                                                                                                   |============================================================================================               |  86%  |                                                                                                                   |=============================================================================================              |  87%  |                                                                                                                   |==============================================================================================             |  88%  |                                                                                                                   |===============================================================================================            |  89%  |                                                                                                                   |================================================================================================           |  90%  |                                                                                                                   |=================================================================================================          |  91%  |                                                                                                                   |==================================================================================================         |  92%  |                                                                                                                   |====================================================================================================       |  93%  |                                                                                                                   |=====================================================================================================      |  94%  |                                                                                                                   |======================================================================================================     |  95%  |                                                                                                                   |=======================================================================================================    |  96%  |                                                                                                                   |========================================================================================================   |  97%  |                                                                                                                   |=========================================================================================================  |  98%  |                                                                                                                   |========================================================================================================== |  99%  |                                                                                                                   |===========================================================================================================| 100%
```

``` r
close(pb)
```

``` r
cat("\nLM:  coverage = ", round(mean(results_lm$coverage), 3), " (SD: ", round(sd(results_lm$coverage), 3), "), width = ", round(mean(results_lm$width), 3), "\n", sep = "")
```

```
## 
## LM:  coverage = 0.907 (SD: 0.058), width = 1.963
```

``` r
cat("GAM: coverage = ", round(mean(results_gam$coverage), 3), " (SD: ", round(sd(results_gam$coverage), 3), "), width = ", round(mean(results_gam$width), 3), "\n", sep = "")
```

```
## GAM: coverage = 0.912 (SD: 0.061), width = 1.848
```

``` r
cat("RF:  coverage = ", round(mean(results_rf$coverage), 3), " (SD: ", round(sd(results_rf$coverage), 3), "), width = ", round(mean(results_rf$width), 3), "\n", sep = "")
```

```
## RF:  coverage = 0.911 (SD: 0.056), width = 2.039
```

``` r
## PART 6: SPATIO-TEMPORAL APPLICATION (NY OZONE DATA)
```

``` r
cat("\n============================================================\n")
```

```
## 
## ============================================================
```

``` r
cat("  PART 6: Spatio-temporal Application (NY Ozone)            \n")
```

```
##   PART 6: Spatio-temporal Application (NY Ozone)
```

``` r
cat("============================================================\n")
```

```
## ============================================================
```

``` r
data("nysptime", package = "bmstdr")
df <- nysptime[complete.cases(nysptime[, c("utmx", "utmy", "y8hrmax", "Day", "Month")]), ]

## Continuous time index (1-62)
df$day_idx <- ifelse(df$Month == 7, df$Day, 31 + df$Day)

s_st <- as.matrix(df[, c("utmx", "utmy")])
y_st <- df$y8hrmax
t_st <- df$day_idx

s_3d <- cbind(s_st, t_st)

pred_fun_gam_3d <- function(s_train, y_train, s_new) {
  train_df <- data.frame(
    x = s_train[, 1], y = s_train[, 2],
    day = s_train[, 3], z = y_train
  )
  fit <- gam(z ~ te(x, y, day, k = c(8, 8, 4)), data = train_df)
  new_df <- data.frame(
    x = s_new[, 1], y = s_new[, 2], day = s_new[, 3]
  )
  as.numeric(predict(fit, newdata = new_df))
}

n_st <- nrow(s_3d)
n_reps_st <- 100
coverages_st <- numeric(n_reps_st)
widths_st    <- numeric(n_reps_st)

pb <- txtProgressBar(min = 0, max = n_reps_st, style = 3)
```

```
##   |                                                                                                                   |                                                                                                           |   0%
```

``` r
for (i in 1:n_reps_st) {
  set.seed(i)
  idx_st <- sample(n_st, floor(0.7 * n_st))
  
  out_st <- scp_geostatistical(
    s_train            = s_3d[idx_st, ],
    y_train            = y_st[idx_st],
    s0                 = s_3d[-idx_st, ],
    pred_fun           = pred_fun_gam_3d,
    t_train            = t_st[idx_st],
    t0                 = t_st[-idx_st],
    temporal_bandwidth = 5,
    alpha              = 0.1,
    split              = 0.5,
    seed               = i
  )
  
  rep_st <- coverage_report(out_st, y_st[-idx_st])
  coverages_st[i] <- rep_st$coverage
  widths_st[i]    <- rep_st$mean_width
  
  setTxtProgressBar(pb, i)
}
```

```
##   |                                                                                                                   |=                                                                                                          |   1%  |                                                                                                                   |==                                                                                                         |   2%  |                                                                                                                   |===                                                                                                        |   3%  |                                                                                                                   |====                                                                                                       |   4%  |                                                                                                                   |=====                                                                                                      |   5%  |                                                                                                                   |======                                                                                                     |   6%  |                                                                                                                   |=======                                                                                                    |   7%  |                                                                                                                   |=========                                                                                                  |   8%  |                                                                                                                   |==========                                                                                                 |   9%  |                                                                                                                   |===========                                                                                                |  10%  |                                                                                                                   |============                                                                                               |  11%  |                                                                                                                   |=============                                                                                              |  12%  |                                                                                                                   |==============                                                                                             |  13%  |                                                                                                                   |===============                                                                                            |  14%  |                                                                                                                   |================                                                                                           |  15%  |                                                                                                                   |=================                                                                                          |  16%  |                                                                                                                   |==================                                                                                         |  17%  |                                                                                                                   |===================                                                                                        |  18%  |                                                                                                                   |====================                                                                                       |  19%  |                                                                                                                   |=====================                                                                                      |  20%  |                                                                                                                   |======================                                                                                     |  21%  |                                                                                                                   |========================                                                                                   |  22%  |                                                                                                                   |=========================                                                                                  |  23%  |                                                                                                                   |==========================                                                                                 |  24%  |                                                                                                                   |===========================                                                                                |  25%  |                                                                                                                   |============================                                                                               |  26%  |                                                                                                                   |=============================                                                                              |  27%  |                                                                                                                   |==============================                                                                             |  28%  |                                                                                                                   |===============================                                                                            |  29%  |                                                                                                                   |================================                                                                           |  30%  |                                                                                                                   |=================================                                                                          |  31%  |                                                                                                                   |==================================                                                                         |  32%  |                                                                                                                   |===================================                                                                        |  33%  |                                                                                                                   |====================================                                                                       |  34%  |                                                                                                                   |=====================================                                                                      |  35%  |                                                                                                                   |=======================================                                                                    |  36%  |                                                                                                                   |========================================                                                                   |  37%  |                                                                                                                   |=========================================                                                                  |  38%  |                                                                                                                   |==========================================                                                                 |  39%  |                                                                                                                   |===========================================                                                                |  40%  |                                                                                                                   |============================================                                                               |  41%  |                                                                                                                   |=============================================                                                              |  42%  |                                                                                                                   |==============================================                                                             |  43%  |                                                                                                                   |===============================================                                                            |  44%  |                                                                                                                   |================================================                                                           |  45%  |                                                                                                                   |=================================================                                                          |  46%  |                                                                                                                   |==================================================                                                         |  47%  |                                                                                                                   |===================================================                                                        |  48%  |                                                                                                                   |====================================================                                                       |  49%  |                                                                                                                   |======================================================                                                     |  50%  |                                                                                                                   |=======================================================                                                    |  51%  |                                                                                                                   |========================================================                                                   |  52%  |                                                                                                                   |=========================================================                                                  |  53%  |                                                                                                                   |==========================================================                                                 |  54%  |                                                                                                                   |===========================================================                                                |  55%  |                                                                                                                   |============================================================                                               |  56%  |                                                                                                                   |=============================================================                                              |  57%  |                                                                                                                   |==============================================================                                             |  58%  |                                                                                                                   |===============================================================                                            |  59%  |                                                                                                                   |================================================================                                           |  60%  |                                                                                                                   |=================================================================                                          |  61%  |                                                                                                                   |==================================================================                                         |  62%  |                                                                                                                   |===================================================================                                        |  63%  |                                                                                                                   |====================================================================                                       |  64%  |                                                                                                                   |======================================================================                                     |  65%  |                                                                                                                   |=======================================================================                                    |  66%  |                                                                                                                   |========================================================================                                   |  67%  |                                                                                                                   |=========================================================================                                  |  68%  |                                                                                                                   |==========================================================================                                 |  69%  |                                                                                                                   |===========================================================================                                |  70%  |                                                                                                                   |============================================================================                               |  71%  |                                                                                                                   |=============================================================================                              |  72%  |                                                                                                                   |==============================================================================                             |  73%  |                                                                                                                   |===============================================================================                            |  74%  |                                                                                                                   |================================================================================                           |  75%  |                                                                                                                   |=================================================================================                          |  76%  |                                                                                                                   |==================================================================================                         |  77%  |                                                                                                                   |===================================================================================                        |  78%  |                                                                                                                   |=====================================================================================                      |  79%  |                                                                                                                   |======================================================================================                     |  80%  |                                                                                                                   |=======================================================================================                    |  81%  |                                                                                                                   |========================================================================================                   |  82%  |                                                                                                                   |=========================================================================================                  |  83%  |                                                                                                                   |==========================================================================================                 |  84%  |                                                                                                                   |===========================================================================================                |  85%  |                                                                                                                   |============================================================================================               |  86%  |                                                                                                                   |=============================================================================================              |  87%  |                                                                                                                   |==============================================================================================             |  88%  |                                                                                                                   |===============================================================================================            |  89%  |                                                                                                                   |================================================================================================           |  90%  |                                                                                                                   |=================================================================================================          |  91%  |                                                                                                                   |==================================================================================================         |  92%  |                                                                                                                   |====================================================================================================       |  93%  |                                                                                                                   |=====================================================================================================      |  94%  |                                                                                                                   |======================================================================================================     |  95%  |                                                                                                                   |=======================================================================================================    |  96%  |                                                                                                                   |========================================================================================================   |  97%  |                                                                                                                   |=========================================================================================================  |  98%  |                                                                                                                   |========================================================================================================== |  99%  |                                                                                                                   |===========================================================================================================| 100%
```

``` r
close(pb)
```

``` r
## Save CSV results
st_results_df <- data.frame(
  replication = 1:n_reps_st,
  coverage    = coverages_st,
  width       = widths_st
)
write.csv(st_results_df, file = "spatio_temporal_results.csv", row.names = FALSE)

cat(sprintf("\nMean coverage (NY Ozone): %.3f (SD: %.3f)\n", mean(coverages_st), sd(coverages_st)))
```

```
## 
## Mean coverage (NY Ozone): 0.897 (SD: 0.017)
```

``` r
cat(sprintf("Mean width (NY Ozone):    %.3f (SD: %.3f)\n", mean(widths_st), sd(widths_st)))
```

```
## Mean width (NY Ozone):    38.273 (SD: 1.293)
```

``` r
## PART 7: FINAL SUMMARY TABLE & SESSION INFO
```

``` r
cat("\n============================================================\n")
```

```
## 
## ============================================================
```

``` r
cat("  FINAL SUMMARY TABLE                                       \n")
```

```
##   FINAL SUMMARY TABLE
```

``` r
cat("============================================================\n")
```

```
## ============================================================
```

``` r
cat("\\begin{table}[htbp]\n")
```

```
## \begin{table}[htbp]
```

``` r
cat("\\centering\n\\small\n")
```

```
## \centering
## \small
```

``` r
cat("\\caption{Summary of empirical results for \\pkg{spconform} on the Meuse and NY ozone datasets.}\n")
```

```
## \caption{Summary of empirical results for \pkg{spconform} on the Meuse and NY ozone datasets.}
```

``` r
cat("\\label{tab:summary_final}\n")
```

```
## \label{tab:summary_final}
```

``` r
cat("\\begin{tabular}{l c c c c}\n\\hline\n")
```

```
## \begin{tabular}{l c c c c}
## \hline
```

``` r
cat("Dataset & Type & $n$ & Nominal coverage & Empirical coverage \\\\\n\\hline\n")
```

```
## Dataset & Type & $n$ & Nominal coverage & Empirical coverage \\
## \hline
```

``` r
cat(sprintf("Meuse (point-referenced) & Geostatistical  & 155        & 0.90  & %.3f \\\\\n", mean(coverages)))
```

```
## Meuse (point-referenced) & Geostatistical  & 155        & 0.90  & 0.920 \\
```

``` r
cat(sprintf("Meuse (aggregated grid)   & Areal           & 21         & 0.80  & %.3f \\\\\n", rep_areal$coverage))
```

```
## Meuse (aggregated grid)   & Areal           & 21         & 0.80  & 0.810 \\
```

``` r
cat(sprintf("NY ozone                 & Spatio-temporal & %d (test) & 0.90* & %.3f \\\\\n", n_st - floor(0.7 * n_st), mean(coverages_st)))
```

```
## NY ozone                 & Spatio-temporal & 514 (test) & 0.90* & 0.897 \\
```

``` r
cat("\\hline\n\\end{tabular}\n\\medskip\n")
```

```
## \hline
## \end{tabular}
## \medskip
```

``` r
cat("\\parbox{\\textwidth}{\\small \\emph{Note:} * The spatio-temporal result is the mean empirical coverage evaluated over 100 independent Monte Carlo data splits ($n_{\\text{test}} = 514$, $\\text{SD} = 0.017$). The mean prediction interval width is $38.273$ ($\\text{SD} = 1.293$).}\n")
```

```
## \parbox{\textwidth}{\small \emph{Note:} * The spatio-temporal result is the mean empirical coverage evaluated over 100 independent Monte Carlo data splits ($n_{\text{test}} = 514$, $\text{SD} = 0.017$). The mean prediction interval width is $38.273$ ($\text{SD} = 1.293$).}
```

``` r
cat("\\end{table}\n\n")
```

```
## \end{table}
```

``` r
cat("Outputs successfully saved:\n")
```

```
## Outputs successfully saved:
```

``` r
cat("  - Figures:      ", file.path(OUTPUT_DIR, "fig1.pdf"), "through", file.path(OUTPUT_DIR, "fig8.pdf"), "\n")
```

```
##   - Figures:       figures/fig1.pdf through figures/fig8.pdf
```

``` r
cat("  - Diagnostics:  spconform_diagnostics.rds\n")
```

```
##   - Diagnostics:  spconform_diagnostics.rds
```

``` r
cat("  - CSV Results:  spatio_temporal_results.csv\n\n")
```

```
##   - CSV Results:  spatio_temporal_results.csv
```

``` r
sessionInfo()
```

```
## R version 4.6.1 (2026-06-24 ucrt)
## Platform: x86_64-w64-mingw32/x64
## Running under: Windows 11 x64 (build 26200)
## 
## Matrix products: default
##   LAPACK version 3.12.1
## 
## locale:
## [1] LC_COLLATE=Arabic_Iraq.utf8  LC_CTYPE=Arabic_Iraq.utf8    LC_MONETARY=Arabic_Iraq.utf8
## [4] LC_NUMERIC=C                 LC_TIME=Arabic_Iraq.utf8    
## 
## time zone: Asia/Baghdad
## tzcode source: internal
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] spconform_0.1.0 bmstdr_0.8.2    Rcpp_1.1.2      ranger_0.18.0   mgcv_1.9-4      nlme_3.1-169    sp_2.2-3       
## 
## loaded via a namespace (and not attached):
##   [1] splines_4.6.1       later_1.4.8         tibble_3.3.1        xts_0.14.2          lifecycle_1.0.5    
##   [6] Rdpack_2.6.6        sf_1.1-2            rstatix_1.1.0       spTimer_3.3.4       rprojroot_2.1.1    
##  [11] StanHeaders_2.32.10 processx_3.9.0      lattice_0.22-9      MASS_7.3-65         crosstalk_1.2.2    
##  [16] backports_1.5.1     magrittr_2.0.5      rmarkdown_2.31      yaml_2.3.12         otel_0.2.0         
##  [21] spam_2.11-4         sessioninfo_1.2.4   pkgbuild_1.4.8      DBI_1.3.0           RColorBrewer_1.1-3 
##  [26] abind_1.4-8         pkgload_1.5.3       purrr_1.2.2         satellite_1.0.6     inline_0.3.21      
##  [31] terra_1.9-34        testthat_3.3.2      units_1.0-1         MatrixModels_0.5-4  pkgdown_2.2.1      
##  [36] codetools_0.2-20    xml2_1.6.0          tidyselect_1.2.1    shape_1.4.6.1       raster_3.6-32      
##  [41] farver_2.1.2        matrixStats_1.5.0   stats4_4.6.1        base64enc_0.1-6     roxygen2_8.1.0     
##  [46] e1071_1.7-17        ellipsis_0.3.3      Formula_1.2-6       survival_3.8-6      iterators_1.0.14   
##  [51] foreach_1.5.2       tools_4.6.1         pak_0.11.1          inlabru_2.15.0      glue_1.8.1         
##  [56] mnormt_2.1.2        gridExtra_2.3.1     xfun_0.60           usethis_3.2.1       dplyr_1.2.1        
##  [61] loo_2.10.1          withr_3.0.3         fastmap_1.2.0       GGally_2.4.0        xopen_1.0.1        
##  [66] boot_1.3-32         SparseM_1.84-2      spData_2.3.5        callr_3.8.0         digest_0.6.39      
##  [71] truncnorm_1.0-9     rcmdcheck_1.4.0     R6_2.6.1            wk_0.9.5            gtools_3.9.5       
##  [76] tidyr_1.3.2         generics_0.1.4      intervals_0.15.5    class_7.3-23        rdtools_0.1.0      
##  [81] prettyunits_1.2.0   htmlwidgets_1.6.4   ggstats_0.13.0      spdep_1.4-2         pkgconfig_2.0.3    
##  [86] gtable_0.3.6        CARBayesdata_3.0    S7_0.2.2            brio_1.1.5          htmltools_0.5.9    
##  [91] carData_3.0-6       spBayes_0.4-9       dotCall64_1.2       fmesher_0.8.0       MCMCpack_1.7-1     
##  [96] scales_1.4.0        png_0.1-9           nanonext_1.10.2     knitr_1.51          rstudioapi_0.19.0  
## [101] coda_0.19-4.1       spacetime_1.3-3     magic_1.6-1         curl_7.1.0          proxy_0.4-29       
## [106] cachem_1.1.0        zoo_1.9-0           KernSmooth_2.23-26  parallel_4.6.1      desc_1.4.3         
## [111] s2_1.1.11           pillar_1.11.1       grid_4.6.1          vctrs_0.7.3         ggpubr_1.0.0       
## [116] mapview_2.11.4      car_3.1-5           CARBayesST_4.0      evaluate_1.0.5      truncdist_1.0-2    
## [121] cli_3.6.6           compiler_4.6.1      rlang_1.3.0         CARBayes_6.1.1      rstantools_2.7.0   
## [126] ggsignif_0.6.4      classInt_0.4-11     ps_1.9.3            fs_2.1.0            rstan_2.32.7       
## [131] deldir_2.0-4        QuickJSR_1.10.0     leaflet_2.2.3       devtools_2.5.2      glmnet_5.0         
## [136] quantreg_6.1        Matrix_1.7-5        leafem_0.2.5        ggplot2_4.0.3       mcmc_0.9-8         
## [141] extraDistr_1.10.0.5 evd_2.3-7.1         rbibutils_2.4.1     igraph_2.3.3        broom_1.0.13       
## [146] memoise_2.0.1       RcppParallel_6.2.0
```

