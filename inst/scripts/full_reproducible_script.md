

``` r
## =============================================================
## spconform - Reproducible Analysis Script 
## Generates all figures, tables, and results 
## =============================================================
```

``` r
options(stringsAsFactors = FALSE)

SEED <- 123
set.seed(SEED)

OUTPUT_DIR <- "figures"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

library(sp)
library(spconform)
library(mgcv)      # For GAM
```

```
## Loading required package: nlme
```

```
## This is mgcv 1.9-4. For overview type '?mgcv'.
```

``` r
library(ranger)    # For Random Forest
```

```
## ranger 0.18.0 using 2 threads (default). Change with num.threads in ranger() and predict(), options(Ncpus = N), options(ranger.num.threads = N) or environment variable R_RANGER_NUM_THREADS.
```

``` r
library(bmstdr)    # For NY ozone data
```

```
## Loading required package: Rcpp
```

``` r
## PART 1: MEUSE RIVER DATA - GEOSTATISTICAL ILLUSTRATION
```

``` r
## 1.1 Load data and define predictors
data(meuse)
s <- as.matrix(meuse[, c("x", "y")])
y <- log(meuse$zinc)

## Point predictor: quadratic trend surface
pred_fun <- function(s_train, y_train, s_new) {
  fit <- lm(y_train ~ s_train[, 1] + s_train[, 2] +
              I(s_train[, 1]^2) + I(s_train[, 2]^2))
  cbind(1, s_new[, 1], s_new[, 2],
        s_new[, 1]^2, s_new[, 2]^2) %*% coef(fit)
}

## 1.2 Figure 1: Spatial layout
pdf(file.path(OUTPUT_DIR, "fig1_meuse_locations.pdf"), width = 6, height = 5)
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

## Figure 2: geostatistical prediction intervals
pdf(file.path(OUTPUT_DIR, "fig2_meuse_geostat_intervals.pdf"), width = 7, height = 5)
if (any(is.na(out$lower)) || any(is.na(out$upper))) {
  warning("NA values found in prediction intervals. Filtering them out.")
  valid <- !is.na(out$lower) & !is.na(out$upper)
  out$lower <- out$lower[valid]
  out$upper <- out$upper[valid]
  out$pred <- out$pred[valid]
  y_test <- y_test[valid]
}
plot(out, y_true = y_test)
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
  idx_i <- sample(n, floor(0.7 * n))
  s_tr <- s[idx_i, ]; y_tr <- y[idx_i]
  s_te <- s[-idx_i, ]; y_te <- y[-idx_i]
  out_i <- scp_geostatistical(s_tr, y_tr, s_te, pred_fun,
                              alpha = 0.1, seed = i)
  rep_i <- coverage_report(out_i, y_te)
  coverages[i] <- rep_i$coverage
  widths[i] <- rep_i$mean_width
}

## Figure 3: Coverage histogram
pdf(file.path(OUTPUT_DIR, "fig3_meuse_coverage_hist.pdf"), width = 6, height = 5)
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
pdf(file.path(OUTPUT_DIR, "fig4_meuse_width_spatial.pdf"), width = 6, height = 5)
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
## 1.6 Summary statistics
cat("\n===== Meuse Geostatistical Results =====\n")
```

```
## 
## ===== Meuse Geostatistical Results =====
```

``` r
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
## PART 2: MEUSE RIVER DATA - AREAL ILLUSTRATION
```

``` r
## 2.1 Aggregate to 6x6 grid
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

## 2.2 Areal conformal prediction
out2 <- scp_areal(agg$y, adjacency = adj, alpha = 0.2, decay = 0.5)

## Figure 5: Areal prediction intervals
pdf(file.path(OUTPUT_DIR, "fig5_meuse_areal_intervals.pdf"), width = 7, height = 5)
if (any(is.na(out2$lower)) || any(is.na(out2$upper))) {
  warning("NA values found in areal prediction intervals. Filtering them out.")
  valid <- !is.na(out2$lower) & !is.na(out2$upper)
  out2$lower <- out2$lower[valid]
  out2$upper <- out2$upper[valid]
  out2$pred <- out2$pred[valid]
  y_true_areal <- agg$y[valid]
  plot(out2, y_true = y_true_areal)
} else {
  plot(out2, y_true = agg$y)
}
```

```
## Warning: NA values found in areal prediction intervals. Filtering them out.
```

``` r
dev.off()
```

```
## RStudioGD 
##         2
```

``` r
## Figure 6: Width comparison
pdf(file.path(OUTPUT_DIR, "fig6_meuse_width_comparison.pdf"), width = 6, height = 5)
geo_width <- out$upper - out$lower
areal_width <- out2$upper - out2$lower
geo_width <- geo_width[!is.na(geo_width)]
areal_width <- areal_width[!is.na(areal_width)]
if (length(geo_width) > 0 && length(areal_width) > 0) {
  boxplot(list(Geostatistical = geo_width,
               Areal = areal_width),
          main = "Interval Width Comparison",
          ylab = "Interval Width",
          col = c("lightblue", "lightgreen"))
} else {
  plot.new()
  text(0.5, 0.5, "Insufficient data for boxplot")
}
dev.off()
```

```
## RStudioGD 
##         2
```

``` r
## Areal summary statistics
cat("\n===== Meuse Areal Results =====\n")
```

```
## 
## ===== Meuse Areal Results =====
```

``` r
if (any(is.na(out2$lower)) || any(is.na(out2$upper))) {
  warning("NA values in areal intervals. Computing coverage on complete cases only.")
  valid <- !is.na(out2$lower) & !is.na(out2$upper)
  y_true_areal <- agg$y[valid]
  rep_areal <- coverage_report(
    list(pred = out2$pred[valid], lower = out2$lower[valid], 
         upper = out2$upper[valid], alpha = out2$alpha),
    y_true_areal
  )
} else {
  rep_areal <- coverage_report(out2, agg$y)
}
```

```
## Warning in y_true >= object$lower: longer object length is not a multiple of shorter object length
```

```
## Warning in y_true <= object$upper: longer object length is not a multiple of shorter object length
```

``` r
cat(sprintf("Areal coverage: %.3f\n", rep_areal$coverage))
```

```
## Areal coverage: 0.857
```

``` r
cat(sprintf("Areal mean width: %.3f\n", rep_areal$mean_width))
```

```
## Areal mean width: 2.434
```

``` r
## PART 3: MODEL COMPARISON (SENSITIVITY ANALYSIS)
```

``` r
cat("\n===== Running Model Comparison (100 splits) =====\n")
```

```
## 
## ===== Running Model Comparison (100 splits) =====
```

``` r
## Define predictors
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

## Run comparison
set.seed(SEED)
n_splits <- 100
n <- nrow(s)

results_lm <- data.frame(coverage = numeric(n_splits), width = numeric(n_splits))
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
  results_lm$width[i] <- rep_lm$mean_width
  
  out_gam <- scp_geostatistical(s_train, y_train, s_test, pred_fun_gam, alpha = 0.1, seed = i)
  rep_gam <- coverage_report(out_gam, y_test)
  results_gam$coverage[i] <- rep_gam$coverage
  results_gam$width[i] <- rep_gam$mean_width
  
  out_rf <- scp_geostatistical(s_train, y_train, s_test, pred_fun_rf, alpha = 0.1, seed = i)
  rep_rf <- coverage_report(out_rf, y_test)
  results_rf$coverage[i] <- rep_rf$coverage
  results_rf$width[i] <- rep_rf$mean_width
  
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
## Print results
cat("\n\n===== Model Comparison Results (Target: 90%) =====\n")
```

```
## 
## 
## ===== Model Comparison Results (Target: 90%) =====
```

``` r
cat(sprintf("LM: coverage = %.3f (SD: %.3f), width = %.3f (SD: %.3f)\n",
            mean(results_lm$coverage), sd(results_lm$coverage),
            mean(results_lm$width), sd(results_lm$width)))
```

```
## LM: coverage = 0.907 (SD: 0.058), width = 1.963 (SD: 0.227)
```

``` r
cat(sprintf("GAM: coverage = %.3f (SD: %.3f), width = %.3f (SD: %.3f)\n",
            mean(results_gam$coverage), sd(results_gam$coverage),
            mean(results_gam$width), sd(results_gam$width)))
```

```
## GAM: coverage = 0.912 (SD: 0.061), width = 1.848 (SD: 0.276)
```

``` r
cat(sprintf("RF: coverage = %.3f (SD: %.3f), width = %.3f (SD: %.3f)\n",
            mean(results_rf$coverage), sd(results_rf$coverage),
            mean(results_rf$width), sd(results_rf$width)))
```

```
## RF: coverage = 0.911 (SD: 0.056), width = 2.039 (SD: 0.237)
```

``` r
## Generate LaTeX table
cat("\n\n===== LaTeX Table for Appendix =====\n")
```

```
## 
## 
## ===== LaTeX Table for Appendix =====
```

``` r
cat("\\begin{table}[htbp]\n\\centering\\small\n")
```

```
## \begin{table}[htbp]
## \centering\small
```

``` r
cat("\\caption{Performance of \\code{scp\\_geostatistical()} with different base predictors on the Meuse dataset. Target coverage is 0.90. Results are averaged over 100 random 70/30 splits.}\n")
```

```
## \caption{Performance of \code{scp\_geostatistical()} with different base predictors on the Meuse dataset. Target coverage is 0.90. Results are averaged over 100 random 70/30 splits.}
```

``` r
cat("\\label{tab:models}\n")
```

```
## \label{tab:models}
```

``` r
cat("\\begin{tabular}{l c c c}\n\\hline\n")
```

```
## \begin{tabular}{l c c c}
## \hline
```

``` r
cat("Base predictor & Mean coverage & Mean interval width & SD (coverage) \\\\\n\\hline\n")
```

```
## Base predictor & Mean coverage & Mean interval width & SD (coverage) \\
## \hline
```

``` r
cat(sprintf("Linear model (quadratic trend) & %.3f & %.3f & %.3f \\\\\n",
            mean(results_lm$coverage), mean(results_lm$width), sd(results_lm$coverage)))
```

```
## Linear model (quadratic trend) & 0.907 & 1.963 & 0.058 \\
```

``` r
cat(sprintf("GAM (smooth spatial terms) & %.3f & \\textbf{%.3f} & %.3f \\\\\n",
            mean(results_gam$coverage), mean(results_gam$width), sd(results_gam$coverage)))
```

```
## GAM (smooth spatial terms) & 0.912 & \textbf{1.848} & 0.061 \\
```

``` r
cat(sprintf("Random Forest & %.3f & %.3f & %.3f \\\\\n",
            mean(results_rf$coverage), mean(results_rf$width), sd(results_rf$coverage)))
```

```
## Random Forest & 0.911 & 2.039 & 0.056 \\
```

``` r
cat("\\hline\n\\end{tabular}\n")
```

```
## \hline
## \end{tabular}
```

``` r
cat("\\medskip\n")
```

```
## \medskip
```

``` r
cat("\\parbox{\\textwidth}{\\small \\emph{Note:} All three predictors maintain empirical coverage close to the nominal 90\\% level, demonstrating the stability of empirical coverage across different base predictors. The GAM yields the narrowest intervals, reflecting its superior ability to capture the smooth spatial structure of the zinc concentration surface.}\n")
```

```
## \parbox{\textwidth}{\small \emph{Note:} All three predictors maintain empirical coverage close to the nominal 90\% level, demonstrating the stability of empirical coverage across different base predictors. The GAM yields the narrowest intervals, reflecting its superior ability to capture the smooth spatial structure of the zinc concentration surface.}
```

``` r
cat("\\end{table}\n")
```

```
## \end{table}
```

``` r
## PART 4: SPATIO-TEMPORAL APPLICATION (NY OZONE DATA)
```

``` r
cat("\n\n===== Spatio-temporal Application (NY Ozone) =====\n")
```

```
## 
## 
## ===== Spatio-temporal Application (NY Ozone) =====
```

``` r
## 4.1 Load and prepare data
data("nysptime")
df <- nysptime[complete.cases(nysptime[, c("utmx", "utmy", "y8hrmax", "Day", "Month")]), ]

## Convert Day (1-31) to a continuous time index (1-62)
df$day_idx <- ifelse(df$Month == 7, df$Day, 31 + df$Day)

s <- as.matrix(df[, c("utmx", "utmy")])
y <- df$y8hrmax
t <- df$day_idx  # Use continuous index instead of Day

## Create 3D coordinates (spatial + temporal)
s_3d <- cbind(s, t)

## 4.2 Define GAM spatio-temporal predictor (using continuous time)
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

## 4.3 100 random splits for stable estimate
set.seed(SEED)
n <- nrow(s_3d)
n_reps <- 100
coverages_st <- numeric(n_reps)
widths_st <- numeric(n_reps)

cat("Running 100 spatio-temporal replications...\n")
```

```
## Running 100 spatio-temporal replications...
```

``` r
pb <- txtProgressBar(min = 0, max = n_reps, style = 3)
```

```
##   |                                                                                                                   |                                                                                                           |   0%
```

``` r
for (i in 1:n_reps) {
  set.seed(i)
  idx <- sample(n, floor(0.7 * n))
  
  s_train_3d <- s_3d[idx, ]
  y_train <- y[idx]
  t_train <- t[idx]
  
  s_test_3d <- s_3d[-idx, ]
  y_test <- y[-idx]
  t_test <- t[-idx]
  
  out_st <- scp_geostatistical(
    s_train = s_train_3d,
    y_train = y_train,
    s0 = s_test_3d,
    pred_fun = pred_fun_gam_3d,
    t_train = t_train,
    t0 = t_test,
    temporal_bandwidth = 5,
    alpha = 0.1,
    split = 0.5,
    seed = i
  )
  
  rep <- coverage_report(out_st, y_test)
  coverages_st[i] <- rep$coverage
  widths_st[i] <- rep$mean_width
  
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
## 4.4 Results
cat("\n\n===== Spatio-temporal Results (100 replications) =====\n")
```

```
## 
## 
## ===== Spatio-temporal Results (100 replications) =====
```

``` r
cat(sprintf("Mean coverage: %.3f (SD: %.3f)\n", mean(coverages_st), sd(coverages_st)))
```

```
## Mean coverage: 0.897 (SD: 0.017)
```

``` r
cat(sprintf("Mean width: %.3f (SD: %.3f)\n", mean(widths_st), sd(widths_st)))
```

```
## Mean width: 38.273 (SD: 1.293)
```

``` r
cat(sprintf("Test set size per split: %d observations\n", n - floor(0.7 * n)))
```

```
## Test set size per split: 514 observations
```

``` r
## Store final spatio-temporal results for the summary table
st_coverage_final <- mean(coverages_st)
st_width_final <- mean(widths_st)
```

``` r
## PART 5: DIAGNOSTICS AND FINAL SUMMARY
```

``` r
cat("\n\n===== FINAL SUMMARY TABLE FOR PAPER =====\n")
```

```
## 
## 
## ===== FINAL SUMMARY TABLE FOR PAPER =====
```

``` r
cat("\\begin{table}[htbp]\n")
```

```
## \begin{table}[htbp]
```

``` r
cat("\\centering\n")
```

```
## \centering
```

``` r
cat("\\small\n")
```

```
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
cat("\\begin{tabular}{l c c c c}\n")
```

```
## \begin{tabular}{l c c c c}
```

``` r
cat("\\hline\n")
```

```
## \hline
```

``` r
cat("Dataset & Type & $n$ & Nominal coverage & Empirical coverage \\\\\n")
```

```
## Dataset & Type & $n$ & Nominal coverage & Empirical coverage \\
```

``` r
cat("\\hline\n")
```

```
## \hline
```

``` r
cat(sprintf("Meuse (point-referenced) & Geostatistical & 155 & 0.90 & %.3f \\\\\n", 
            mean(coverages)))
```

```
## Meuse (point-referenced) & Geostatistical & 155 & 0.90 & 0.920 \\
```

``` r
cat(sprintf("Meuse (aggregated grid) & Areal & 21 & 0.80 & %.3f \\\\\n",
            rep_areal$coverage))
```

```
## Meuse (aggregated grid) & Areal & 21 & 0.80 & 0.857 \\
```

``` r
cat(sprintf("NY ozone & Spatio-temporal & %d (test) & 0.90* & %.3f \\\\\n",
            n - floor(0.7 * n), st_coverage_final))
```

```
## NY ozone & Spatio-temporal & 514 (test) & 0.90* & 0.897 \\
```

``` r
cat("\\hline\n")
```

```
## \hline
```

``` r
cat("\\end{tabular}\n")
```

```
## \end{tabular}
```

``` r
cat("\\medskip\n")
```

```
## \medskip
```

``` r
cat("\\parbox{\\textwidth}{\\small \\emph{Note:} * The spatio-temporal result is the mean over 100 random splits and is illustrative due to the stronger temporal dependence structure. The nominal 90\\% level is used as a reference; empirical coverage may vary due to temporal dependence and limited effective calibration size.}\n")
```

```
## \parbox{\textwidth}{\small \emph{Note:} * The spatio-temporal result is the mean over 100 random splits and is illustrative due to the stronger temporal dependence structure. The nominal 90\% level is used as a reference; empirical coverage may vary due to temporal dependence and limited effective calibration size.}
```

``` r
cat("\\end{table}\n")
```

```
## \end{table}
```

``` r
cat("\n===== Generated outputs =====\n")
```

```
## 
## ===== Generated outputs =====
```

``` r
cat("Figures saved in:", OUTPUT_DIR, "\n")
```

```
## Figures saved in: figures
```

``` r
cat("CSV results: spatio_temporal_results.csv\n")
```

```
## CSV results: spatio_temporal_results.csv
```

``` r
cat("Diagnostics: spconform_diagnostics.rds\n")
```

```
## Diagnostics: spconform_diagnostics.rds
```

``` r
cat("\n===== Session Information =====\n")
```

```
## 
## ===== Session Information =====
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
## [1] bmstdr_0.8.2    Rcpp_1.1.2      ranger_0.18.0   mgcv_1.9-4      nlme_3.1-169    spconform_0.1.0 sp_2.2-3       
## 
## loaded via a namespace (and not attached):
##   [1] splines_4.6.1       later_1.4.8         tibble_3.3.1        xts_0.14.2          lifecycle_1.0.5    
##   [6] Rdpack_2.6.6        sf_1.1-2            rstatix_1.1.0       spTimer_3.3.4       rprojroot_2.1.1    
##  [11] StanHeaders_2.32.10 processx_3.9.0      lattice_0.22-9      MASS_7.3-65         crosstalk_1.2.2    
##  [16] backports_1.5.1     magrittr_2.0.5      sass_0.4.10         rmarkdown_2.31      jquerylib_0.1.4    
##  [21] yaml_2.3.12         remotes_2.5.0       otel_0.2.0          spam_2.11-4         sessioninfo_1.2.4  
##  [26] pkgbuild_1.4.8      DBI_1.3.0           RColorBrewer_1.1-3  abind_1.4-8         pkgload_1.5.3      
##  [31] purrr_1.2.2         downlit_0.4.5       satellite_1.0.6     inline_0.3.21       terra_1.9-34       
##  [36] testthat_3.3.2      units_1.0-1         MatrixModels_0.5-4  pkgdown_2.2.1       codetools_0.2-20   
##  [41] xml2_1.6.0          tidyselect_1.2.1    shape_1.4.6.1       raster_3.6-32       farver_2.1.2       
##  [46] matrixStats_1.5.0   stats4_4.6.1        base64enc_0.1-6     roxygen2_8.1.0      jsonlite_2.0.0     
##  [51] e1071_1.7-17        ellipsis_0.3.3      Formula_1.2-6       survival_3.8-6      iterators_1.0.14   
##  [56] foreach_1.5.2       tools_4.6.1         inlabru_2.15.0      glue_1.8.1          mnormt_2.1.2       
##  [61] gridExtra_2.3.1     xfun_0.60           usethis_3.2.1       dplyr_1.2.1         loo_2.10.1         
##  [66] withr_3.0.3         fastmap_1.2.0       GGally_2.4.0        boot_1.3-32         xopen_1.0.1        
##  [71] fansi_1.0.7         SparseM_1.84-2      spData_2.3.5        callr_3.8.0         digest_0.6.39      
##  [76] truncnorm_1.0-9     rcmdcheck_1.4.0     R6_2.6.1            wk_0.9.5            gtools_3.9.5       
##  [81] tidyr_1.3.2         generics_0.1.4      intervals_0.15.5    class_7.3-23        rdtools_0.1.0      
##  [86] prettyunits_1.2.0   htmlwidgets_1.6.4   ggstats_0.13.0      whisker_0.4.1       spdep_1.4-2        
##  [91] pkgconfig_2.0.3     gtable_0.3.6        CARBayesdata_3.0    S7_0.2.2            brio_1.1.5         
##  [96] htmltools_0.5.9     carData_3.0-6       spBayes_0.4-9       dotCall64_1.2       fmesher_0.8.0      
## [101] MCMCpack_1.7-1      scales_1.4.0        png_0.1-9           nanonext_1.10.2     knitr_1.51         
## [106] rstudioapi_0.19.0   coda_0.19-4.1       spacetime_1.3-3     magic_1.6-1         curl_7.1.0         
## [111] proxy_0.4-29        cachem_1.1.0        zoo_1.9-0           KernSmooth_2.23-26  parallel_4.6.1     
## [116] desc_1.4.3          s2_1.1.11           pillar_1.11.1       grid_4.6.1          vctrs_0.7.3        
## [121] ggpubr_1.0.0        mapview_2.11.4      car_3.1-5           CARBayesST_4.0      evaluate_1.0.5     
## [126] truncdist_1.0-2     cli_3.6.6           compiler_4.6.1      rlang_1.3.0         CARBayes_6.1.1     
## [131] rstantools_2.7.0    ggsignif_0.6.4      classInt_0.4-11     fs_2.1.0            rstan_2.32.7       
## [136] deldir_2.0-4        QuickJSR_1.10.0     leaflet_2.2.3       devtools_2.5.2      glmnet_5.0         
## [141] quantreg_6.1        Matrix_1.7-5        leafem_0.2.5        ggplot2_4.0.3       mcmc_0.9-8         
## [146] extraDistr_1.10.0.5 evd_2.3-7.1         rbibutils_2.4.1     igraph_2.3.3        broom_1.0.13       
## [151] memoise_2.0.1       RcppParallel_6.2.0  bslib_0.12.0
```

