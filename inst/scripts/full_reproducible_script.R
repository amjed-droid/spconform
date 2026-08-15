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

## -----------------------------------------------------------------------------
## PART 1: MEUSE RIVER DATA - GEOSTATISTICAL ILLUSTRATION
## -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("  PART 1: Meuse River Data - Geostatistical Analysis        \n")
cat("============================================================\n")

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

## 1.5 Figure 4: Spatial distribution of interval width
plot_df <- data.frame(x = s_test[, 1], y = s_test[, 2],
                      width = out$upper - out$lower)
pdf(file.path(OUTPUT_DIR, "fig4.pdf"), width = 6, height = 5)
plot(plot_df$x, plot_df$y, cex = plot_df$width, pch = 19,
     col = rgb(0.2, 0.4, 0.8, 0.5),
     xlab = "X coordinate", ylab = "Y coordinate",
     main = "Spatial Distribution of Interval Width")
dev.off()

## 1.6 Summary statistics for Geostatistical section
cat(sprintf("Single split coverage: %.3f\n", coverage_report(out, y_test)$coverage))
cat(sprintf("Mean coverage (50 splits): %.3f (SD: %.3f)\n", 
            mean(coverages), sd(coverages)))
cat(sprintf("Mean width (50 splits): %.3f\n", mean(widths)))


## -----------------------------------------------------------------------------
## PART 2: DIAGNOSTIC REPORT (FIGURE 5)
## -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("  PART 2: Spatial Diagnostics (Figure 5)                    \n")
cat("============================================================\n")

pdf(file.path(OUTPUT_DIR, "fig8.pdf"), width = 8.5, height = 7)
diag_meuse <- diagnose(
  object  = out,
  y_true  = y_test,
  s_test  = s_test,
  n_bins  = 4,
  plot    = TRUE
)
dev.off()

saveRDS(diag_meuse, file = "spconform_diagnostics.rds")
cat("Figure 5 (diagnostics) saved to:", file.path(OUTPUT_DIR, "fig8.pdf"), "\n")


## -----------------------------------------------------------------------------
## PART 3: MEUSE RIVER DATA - AREAL ILLUSTRATION
## -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("  PART 3: Meuse River Data - Areal Lattice Analysis         \n")
cat("============================================================\n")

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
cat(sprintf("Areal mean width: %.3f\n", rep_areal$mean_width))


## -----------------------------------------------------------------------------
## PART 4: HELD-OUT AREAL EVALUATION (FIGURE 8)
## -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("  PART 4: Held-out Areal Evaluation (Figure 8)              \n")
cat("============================================================\n")

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


## -----------------------------------------------------------------------------
## PART 5: SENSITIVITY ANALYSIS: MODEL COMPARISON (100 SPLITS)
## -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("  PART 5: Sensitivity Analysis (Linear vs GAM vs RF)        \n")
cat("============================================================\n")

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
close(pb)

cat("\nLM:  coverage = ", round(mean(results_lm$coverage), 3), " (SD: ", round(sd(results_lm$coverage), 3), "), width = ", round(mean(results_lm$width), 3), "\n", sep = "")
cat("GAM: coverage = ", round(mean(results_gam$coverage), 3), " (SD: ", round(sd(results_gam$coverage), 3), "), width = ", round(mean(results_gam$width), 3), "\n", sep = "")
cat("RF:  coverage = ", round(mean(results_rf$coverage), 3), " (SD: ", round(sd(results_rf$coverage), 3), "), width = ", round(mean(results_rf$width), 3), "\n", sep = "")


## -----------------------------------------------------------------------------
## PART 6: SPATIO-TEMPORAL APPLICATION (NY OZONE DATA)
## -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("  PART 6: Spatio-temporal Application (NY Ozone)            \n")
cat("============================================================\n")

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
close(pb)

## Save CSV results
st_results_df <- data.frame(
  replication = 1:n_reps_st,
  coverage    = coverages_st,
  width       = widths_st
)
write.csv(st_results_df, file = "spatio_temporal_results.csv", row.names = FALSE)

cat(sprintf("\nMean coverage (NY Ozone): %.3f (SD: %.3f)\n", mean(coverages_st), sd(coverages_st)))
cat(sprintf("Mean width (NY Ozone):    %.3f (SD: %.3f)\n", mean(widths_st), sd(widths_st)))


## -----------------------------------------------------------------------------
## PART 7: FINAL SUMMARY TABLE & SESSION INFO
## -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("  FINAL SUMMARY TABLE                                       \n")
cat("============================================================\n")

cat("\\begin{table}[htbp]\n")
cat("\\centering\n\\small\n")
cat("\\caption{Summary of empirical results for \\pkg{spconform} on the Meuse and NY ozone datasets.}\n")
cat("\\label{tab:summary_final}\n")
cat("\\begin{tabular}{l c c c c}\n\\hline\n")
cat("Dataset & Type & $n$ & Nominal coverage & Empirical coverage \\\\\n\\hline\n")
cat(sprintf("Meuse (point-referenced) & Geostatistical  & 155        & 0.90  & %.3f \\\\\n", mean(coverages)))
cat(sprintf("Meuse (aggregated grid)   & Areal           & 21         & 0.80  & %.3f \\\\\n", rep_areal$coverage))
cat(sprintf("NY ozone                 & Spatio-temporal & %d (test) & 0.90* & %.3f \\\\\n", n_st - floor(0.7 * n_st), mean(coverages_st)))
cat("\\hline\n\\end{tabular}\n\\medskip\n")
cat("\\parbox{\\textwidth}{\\small \\emph{Note:} * The spatio-temporal result is the mean empirical coverage evaluated over 100 independent Monte Carlo data splits ($n_{\\text{test}} = 514$, $\\text{SD} = 0.017$). The mean prediction interval width is $38.273$ ($\\text{SD} = 1.293$).}\n")
cat("\\end{table}\n\n")

cat("Outputs successfully saved:\n")
cat("  - Figures:      ", file.path(OUTPUT_DIR, "fig1.pdf"), "through", file.path(OUTPUT_DIR, "fig8.pdf"), "\n")
cat("  - Diagnostics:  spconform_diagnostics.rds\n")
cat("  - CSV Results:  spatio_temporal_results.csv\n\n")

sessionInfo()
