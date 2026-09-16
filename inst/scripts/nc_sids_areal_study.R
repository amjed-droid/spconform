## =============================================================================
## Script: nc_sids_areal_study.R
## Title: Large-Scale Real Areal Benchmark: North Carolina SIDS (n = 100)
## Description: Standalone, publication-ready areal conformal prediction study
##              using the canonical North Carolina Sudden Infant Death Syndrome
##              (SIDS) benchmark dataset (n=100 counties), replacing small-lattice
##              toy examples (e.g. 6x6 grid, n=21) for JSS submission.
##
## Package: spconform
## Author: Ahmed Sattar Jabbar
## License: GPL (>= 3)
## =============================================================================

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(spconform)
  }
  library(stats)
  library(graphics)
  library(grDevices)
})

# Verify required spatial packages for extracting the benchmark shapefile
if (!requireNamespace("sf", quietly = TRUE) || !requireNamespace("spdep", quietly = TRUE)) {
  stop("This script requires packages 'sf' and 'spdep'. Please install them via install.packages(c('sf', 'spdep')).")
}

library(sf)
library(spdep)

cat("\n===============================================================================\n")
cat("  SPCONFORM: Areal Benchmark Analysis on North Carolina SIDS (n = 100)        \n")
cat("  Prepared for Journal of Statistical Software (JSS) Manuscript              \n")
cat("===============================================================================\n\n")

## Output directory setup (CRAN-compliant default to tempdir, configurable via env)
OUTPUT_DIR <- Sys.getenv("SPCONFORM_OUTPUT_DIR", unset = file.path(tempdir(), "nc_sids_outputs"))
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
cat(sprintf("[Setup] Output directory for figures and tables: %s\n\n", OUTPUT_DIR))

## -----------------------------------------------------------------------------
## 1. LOAD DATA & CONSTRUCT REAL GRAPH TOPOLOGY
## -----------------------------------------------------------------------------
cat("--- Step 1: Loading North Carolina SIDS Dataset (n = 100) ---\n")
nc_path <- system.file("shape/nc.shp", package = "sf")
if (!file.exists(nc_path)) {
  stop("North Carolina shapefile not found in sf package installation.")
}

nc <- sf::st_read(nc_path, quiet = TRUE)
n_counties <- nrow(nc)
cat(sprintf("Loaded %d counties from North Carolina.\n", n_counties))

## Centroids for spatial coordinates
suppressWarnings({
  coords <- sf::st_coordinates(sf::st_centroid(nc))
})
colnames(coords) <- c("long", "lat")

## Build Queen contiguity graph
nb_queen <- spdep::poly2nb(nc, queen = TRUE)
adj_mat  <- spdep::nb2mat(nb_queen, style = "B", zero.policy = TRUE)
colnames(adj_mat) <- rownames(adj_mat) <- nc$NAME

n_links <- sum(adj_mat) / 2
degrees <- colSums(adj_mat)
cat(sprintf("Adjacency graph constructed: %d nodes, %d undirected edges (mean degree: %.2f, range: %d - %d).\n",
            n_counties, n_links, mean(degrees), min(degrees), max(degrees)))

## Define epidemiologic response: log-transformed SIDS rate (1974-1978)
## Log rate per 1,000 live births (with continuity correction +1)
y_sids <- log(1000 * (nc$SID74 + 1) / nc$BIR74)
cov_nwprop <- matrix(nc$NWBIR74 / nc$BIR74, ncol = 1)  # Non-white birth proportion (known spatial covariate)

cat(sprintf("Response variable: log SIDS rate (mean = %.3f, sd = %.3f, min = %.3f, max = %.3f)\n\n",
            mean(y_sids), sd(y_sids), min(y_sids), max(y_sids)))

## -----------------------------------------------------------------------------
## 2. FULL CONFORMAL CALIBRATION ACROSS MULTIPLE COVERAGE TARGETS
## -----------------------------------------------------------------------------
cat("--- Step 2: Evaluating scp_areal across nominal coverage targets ---\n")

alphas <- c(0.20, 0.10, 0.05)
coverage_targets <- c("80%", "90%", "95%")

eval_summary <- data.frame(
  Target_Coverage = coverage_targets,
  Alpha           = alphas,
  Empirical_Cov   = numeric(length(alphas)),
  Covered_Count   = integer(length(alphas)),
  Total_Count     = rep(n_counties, length(alphas)),
  Mean_Width      = numeric(length(alphas)),
  Median_Width    = numeric(length(alphas)),
  IQR_Width       = numeric(length(alphas)),
  Min_Width       = numeric(length(alphas)),
  Max_Width       = numeric(length(alphas))
)

fitted_models <- list()

for (k in seq_along(alphas)) {
  a <- alphas[k]
  fit <- scp_areal(y = y_sids, adjacency = adj_mat, alpha = a, decay = 0.5)
  fitted_models[[as.character(a)]] <- fit
  
  covered <- (y_sids >= fit$lower) & (y_sids <= fit$upper)
  widths  <- fit$upper - fit$lower
  
  eval_summary$Empirical_Cov[k] <- mean(covered, na.rm = TRUE)
  eval_summary$Covered_Count[k] <- sum(covered, na.rm = TRUE)
  eval_summary$Mean_Width[k]    <- mean(widths, na.rm = TRUE)
  eval_summary$Median_Width[k]  <- median(widths, na.rm = TRUE)
  eval_summary$IQR_Width[k]     <- IQR(widths, na.rm = TRUE)
  eval_summary$Min_Width[k]     <- min(widths, na.rm = TRUE)
  eval_summary$Max_Width[k]     <- max(widths, na.rm = TRUE)
}

print(eval_summary, row.names = FALSE)
cat("\n")

## -----------------------------------------------------------------------------
## 3. COMPARISON WITH COVARIATE-AUGMENTED AREAL PREDICTOR
## -----------------------------------------------------------------------------
cat("--- Step 3: Comparing Default Predictor vs Covariate-Augmented Predictor ---\n")

## Custom predictor incorporating non-white birth proportion covariate
pred_covariate <- function(y_train, X_train, idx_train, idx_target, adjacency) {
  nw_train  <- as.numeric(X_train[, 1])
  nw_target <- as.numeric(cov_nwprop[idx_target, 1])
  fit_lm    <- lm(y_train ~ nw_train)
  as.numeric(predict(fit_lm, newdata = data.frame(nw_train = nw_target)))
}

fit_cov <- scp_areal(y = y_sids, X = cov_nwprop, adjacency = adj_mat,
                     pred_fun = pred_covariate, alpha = 0.10, decay = 0.5)

cov_covered <- (y_sids >= fit_cov$lower) & (y_sids <= fit_cov$upper)
cov_widths  <- fit_cov$upper - fit_cov$lower

cat(sprintf("Default Predictor (Neighborhood Mean, alpha=0.10):       Coverage = %.2f%%, Mean Width = %.3f\n",
            eval_summary$Empirical_Cov[eval_summary$Alpha == 0.10] * 100,
            eval_summary$Mean_Width[eval_summary$Alpha == 0.10]))
cat(sprintf("Covariate-Augmented Predictor (LM with NWBIR, alpha=0.10): Coverage = %.2f%%, Mean Width = %.3f\n\n",
            mean(cov_covered) * 100, mean(cov_widths)))

## -----------------------------------------------------------------------------
## 4. SPATIAL STRATIFICATION & DEGREE SENSITIVITY ANALYSIS
## -----------------------------------------------------------------------------
cat("--- Step 4: Spatial Stratification and Topological Connectivity Analysis ---\n")

fit_90 <- fitted_models[["0.1"]]
covered_90 <- (y_sids >= fit_90$lower) & (y_sids <= fit_90$upper)
widths_90  <- fit_90$upper - fit_90$lower

# Degree stratification (Sparse <= 3 neighbors, Moderate 4-6, Dense >= 7)
degree_cat <- cut(degrees, breaks = c(0, 3, 6, 20), labels = c("Sparse (<=3)", "Moderate (4-6)", "Dense (>=7)"))
degree_table <- aggregate(data.frame(Coverage = covered_90, Width = widths_90),
                          by = list(Topology = degree_cat), FUN = mean)
degree_table$Count <- as.vector(table(degree_cat))
cat("Empirical Performance by Topological Node Degree:\n")
print(degree_table)
cat("\n")

# Comprehensive spatial diagnostic report using centroids
cat("Running full spatial diagnostics via diagnose()...\n")
diag_report <- diagnose(fit_90, y_true = y_sids, s_test = coords, n_bins = 4, plot = FALSE)
cat(sprintf("Marginal Coverage: %.2f (Nominal: 0.90)\n", diag_report$marginal$coverage))
cat(sprintf("Boundary Effect: Near boundary coverage = %.4f | Core coverage = %.4f\n\n",
            diag_report$boundary$near_boundary$coverage,
            diag_report$boundary$far_boundary$coverage))

## -----------------------------------------------------------------------------
## 5. MONTE CARLO CROSS-VALIDATION (100 REPETITIONS, 70/30 TRAIN-TEST SPLITS)
## -----------------------------------------------------------------------------
cat("--- Step 5: Monte Carlo Cross-Validation (100 Out-of-Sample Repetitions) ---\n")

set.seed(42)
n_mc_reps <- 100
mc_results <- data.frame(
  Rep              = seq_len(n_mc_reps),
  Train_Coverage   = numeric(n_mc_reps),
  Train_Width      = numeric(n_mc_reps),
  Test_Coverage    = numeric(n_mc_reps),
  Test_Width       = numeric(n_mc_reps)
)

for (r in seq_len(n_mc_reps)) {
  # 70% train (70 counties), 30% test (30 counties)
  train_idx <- sample(n_counties, size = 70)
  test_idx  <- setdiff(seq_len(n_counties), train_idx)
  
  y_train <- y_sids[train_idx]
  y_test  <- y_sids[test_idx]
  adj_train <- adj_mat[train_idx, train_idx]
  
  # Train LOO conformal model
  fit_tr <- suppressWarnings(
    scp_areal(y = y_train, adjacency = adj_train, alpha = 0.10, decay = 0.5)
  )
  tr_cov <- (y_train >= fit_tr$lower) & (y_train <= fit_tr$upper)
  mc_results$Train_Coverage[r] <- mean(tr_cov, na.rm = TRUE)
  mc_results$Train_Width[r]    <- mean(fit_tr$upper - fit_tr$lower, na.rm = TRUE)
  
  # Predict on unseen test counties using graph-weighted conformal calibration
  scores_train <- abs(y_train - fit_tr$pred)
  
  test_cov_vec <- logical(length(test_idx))
  test_width_vec <- numeric(length(test_idx))
  
  for (j in seq_along(test_idx)) {
    te_id <- test_idx[j]
    
    # Graph weights from test unit to training units
    w_all <- areal_neighbor_weights(te_id, adj_mat, decay = 0.5)
    w_tr  <- w_all[train_idx]
    
    if (sum(w_tr) > 0) {
      w_tr_norm <- w_tr / sum(w_tr)
      pred_te   <- sum(w_tr_norm * y_train)
    } else {
      w_tr_norm <- rep(1 / length(y_train), length(y_train))
      pred_te   <- mean(y_train)
    }
    
    # Conformal quantile at finite-sample level tau = ceil((1-alpha)*(m+1))/m
    m_tr  <- length(scores_train)
    tau   <- min(1, ceiling((1 - 0.10) * (m_tr + 1)) / m_tr)
    ord   <- order(scores_train)
    cum_w <- cumsum(w_tr_norm[ord])
    q_idx <- which(cum_w >= tau)[1]
    if (is.na(q_idx)) q_idx <- length(scores_train)
    q_val <- scores_train[ord][q_idx]
    
    lower_te <- pred_te - q_val
    upper_te <- pred_te + q_val
    
    test_cov_vec[j] <- (y_test[j] >= lower_te) && (y_test[j] <= upper_te)
    test_width_vec[j] <- upper_te - lower_te
  }
  
  mc_results$Test_Coverage[r] <- mean(test_cov_vec)
  mc_results$Test_Width[r]    <- mean(test_width_vec)
}

cat(sprintf("Monte Carlo Results across %d splits (Nominal Target: 90%%):\n", n_mc_reps))
cat(sprintf("  Mean Test Coverage:   %.2f%% (SD: %.2f%%, 95%% CI: [%.2f%%, %.2f%%])\n",
            mean(mc_results$Test_Coverage) * 100, sd(mc_results$Test_Coverage) * 100,
            (mean(mc_results$Test_Coverage) - 1.96 * sd(mc_results$Test_Coverage) / sqrt(n_mc_reps)) * 100,
            (mean(mc_results$Test_Coverage) + 1.96 * sd(mc_results$Test_Coverage) / sqrt(n_mc_reps)) * 100))
cat(sprintf("  Mean Test Width:      %.3f (SD: %.3f)\n\n",
            mean(mc_results$Test_Width), sd(mc_results$Test_Width)))

## -----------------------------------------------------------------------------
## 6. GENERATE PUBLICATION-QUALITY FIGURES FOR JSS
## -----------------------------------------------------------------------------
cat("--- Step 6: Generating Publication-Ready Figures (PDF & PNG) ---\n")

## Figure 1: Prediction intervals for all 100 counties
for (ext in c("pdf", "png")) {
  fig_file <- file.path(OUTPUT_DIR, paste0("fig_nc_intervals.", ext))
  if (ext == "pdf") pdf(fig_file, width = 9, height = 5)
  if (ext == "png") png(fig_file, width = 1800, height = 1000, res = 200)
  
  ord <- order(y_sids)
  par(mar = c(4.2, 4.2, 2.5, 1.2))
  plot(seq_len(n_counties), y_sids[ord], type = "n",
       ylim = range(c(fit_90$lower, fit_90$upper)),
       xlab = "North Carolina Counties (Ordered by Observed SIDS Rate)",
       ylab = "Log SIDS Rate (per 1,000 live births)",
       main = "Localized Conformal Prediction Intervals on NC SIDS (n = 100, Nominal 90%)")
  
  grid(col = "gray88", lty = "dotted")
  
  for (i in seq_len(n_counties)) {
    orig_i <- ord[i]
    is_cov <- covered_90[orig_i]
    col_line <- if (is_cov) "#2B8CBE" else "#E41A1C"
    lwd_line <- if (is_cov) 1.2 else 2.2
    
    segments(i, fit_90$lower[orig_i], i, fit_90$upper[orig_i], col = col_line, lwd = lwd_line)
    points(i, fit_90$pred[orig_i], pch = 3, cex = 0.5, col = "gray40")
  }
  points(seq_len(n_counties), y_sids[ord], pch = 19, cex = 0.6, col = "black")
  
  legend("topleft",
         legend = c(sprintf("Covered Interval (n=%d)", sum(covered_90)),
                    sprintf("Uncovered Interval (n=%d)", sum(!covered_90)),
                    "Point Prediction", "Observed Rate"),
         col = c("#2B8CBE", "#E41A1C", "gray40", "black"),
         lwd = c(1.5, 2.2, NA, NA),
         pch = c(NA, NA, 3, 19),
         pt.cex = c(NA, NA, 0.7, 0.8),
         bty = "o", bg = "white")
  dev.off()
}

## Figure 2: Spatial Map of Prediction Intervals & Coverage
for (ext in c("pdf", "png")) {
  fig_file <- file.path(OUTPUT_DIR, paste0("fig_nc_spatial_map.", ext))
  if (ext == "pdf") pdf(fig_file, width = 8.5, height = 4.5)
  if (ext == "png") png(fig_file, width = 1800, height = 950, res = 200)
  
  par(mar = c(1, 1, 2.5, 1))
  nc$covered_status <- factor(covered_90, levels = c(TRUE, FALSE), labels = c("Covered (90%)", "Uncovered"))
  
  # Color palette by interval width
  n_shades <- 5
  pal <- hcl.colors(n_shades, palette = "Viridis", rev = TRUE)
  width_bins <- cut(widths_90, breaks = quantile(widths_90, probs = seq(0, 1, length.out = n_shades + 1)),
                    include.lowest = TRUE)
  
  plot(sf::st_geometry(nc), col = pal[width_bins], border = "gray30", lwd = 0.6,
       main = "North Carolina SIDS: Conformal Interval Widths & Coverage (Nominal 90%)")
  
  # Highlight uncovered counties in red border
  uncovered_ids <- which(!covered_90)
  if (length(uncovered_ids) > 0) {
    plot(sf::st_geometry(nc[uncovered_ids, ]), col = "red", border = "darkred", lwd = 1.8, add = TRUE)
  }
  
  legend("bottomleft",
         legend = c("Uncovered County", levels(width_bins)),
         fill = c("red", pal),
         border = "gray30",
         title = "Interval Width / Status",
         bty = "o", bg = "white", cex = 0.8)
  dev.off()
}

## Figure 3: Interval Width vs Node Degree
for (ext in c("pdf", "png")) {
  fig_file <- file.path(OUTPUT_DIR, paste0("fig_nc_width_vs_degree.", ext))
  if (ext == "pdf") pdf(fig_file, width = 6.5, height = 5)
  if (ext == "png") png(fig_file, width = 1400, height = 1100, res = 200)
  
  par(mar = c(4.2, 4.2, 2.5, 1.2))
  boxplot(widths_90 ~ degrees,
          xlab = "County Degree (Number of Neighboring Counties)",
          ylab = "Conformal Prediction Interval Width",
          main = "Uncertainty Adaptation: Interval Width vs Graph Degree",
          col = "#A6BDDB", border = "#02818A")
  grid(col = "gray88", lty = "dotted")
  dev.off()
}

## Figure 4: Monte Carlo Out-of-Sample Coverage Distribution
for (ext in c("pdf", "png")) {
  fig_file <- file.path(OUTPUT_DIR, paste0("fig_nc_monte_carlo.", ext))
  if (ext == "pdf") pdf(fig_file, width = 7, height = 5)
  if (ext == "png") png(fig_file, width = 1500, height = 1100, res = 200)
  
  par(mar = c(4.2, 4.2, 2.5, 1.2))
  hist(mc_results$Test_Coverage * 100, breaks = 15, col = "#74C476", border = "#238B45",
       xlab = "Empirical Test Coverage (%)",
       ylab = "Number of Splits (out of 100)",
       main = "Out-of-Sample Conformal Coverage Across 100 Random Splits")
  abline(v = 90, col = "red", lwd = 2.5, lty = 2)
  abline(v = mean(mc_results$Test_Coverage) * 100, col = "darkblue", lwd = 2)
  legend("topright",
         legend = c("Nominal Target (90%)",
                    sprintf("Mean Empirical (%.1f%%)", mean(mc_results$Test_Coverage) * 100)),
         col = c("red", "darkblue"), lty = c(2, 1), lwd = c(2.5, 2), bty = "o", bg = "white")
  dev.off()
}

cat("All publication figures successfully saved.\n\n")

## -----------------------------------------------------------------------------
## 7. SAVE DETAILED DATA TABLES (CSV) FOR JSS PAPER
## -----------------------------------------------------------------------------
cat("--- Step 7: Exporting Data Tables for JSS Manuscript ---\n")

## County-level details table
county_df <- data.frame(
  County     = nc$NAME,
  Observed_Y = round(y_sids, 4),
  Predicted  = round(fit_90$pred, 4),
  Lower_90   = round(fit_90$lower, 4),
  Upper_90   = round(fit_90$upper, 4),
  Width_90   = round(widths_90, 4),
  Covered_90 = covered_90,
  Degree     = degrees
)
write.csv(county_df, file.path(OUTPUT_DIR, "nc_sids_county_intervals.csv"), row.names = FALSE)

## Evaluation summary table
write.csv(eval_summary, file.path(OUTPUT_DIR, "nc_sids_coverage_summary.csv"), row.names = FALSE)

## Monte Carlo 100 splits table
write.csv(mc_results, file.path(OUTPUT_DIR, "nc_sids_monte_carlo_results.csv"), row.names = FALSE)

cat(sprintf("Saved tables to: %s\n", OUTPUT_DIR))
cat("  - nc_sids_county_intervals.csv\n")
cat("  - nc_sids_coverage_summary.csv\n")
cat("  - nc_sids_monte_carlo_results.csv\n\n")

## -----------------------------------------------------------------------------
## 8. JSS LATEX/MARKDOWN TABLE SNIPPET (Ready to copy-paste)
## -----------------------------------------------------------------------------
cat("===============================================================================\n")
cat("  MARKDOWN TABLE READY FOR JSS MANUSCRIPT / RESPONSE TO REVIEWERS             \n")
cat("===============================================================================\n\n")

cat("| Benchmark Dataset | Type | Units (n) | Target Coverage | Empirical Coverage | Mean Interval Width | Median Width |\n")
cat("|:---|:---|:---:|:---:|:---:|:---:|:---:|\n")
for (k in seq_along(alphas)) {
  cat(sprintf("| North Carolina SIDS | Areal (Counties) | %d | %s | %.2f%% (%d/%d) | %.3f | %.3f |\n",
              n_counties, eval_summary$Target_Coverage[k],
              eval_summary$Empirical_Cov[k] * 100,
              eval_summary$Covered_Count[k], n_counties,
              eval_summary$Mean_Width[k], eval_summary$Median_Width[k]))
}
cat(sprintf("| Monte Carlo (100 splits) | Held-Out Areal | %d train / %d test | 90%% | %.2f%% (± %.2f%%) | %.3f | %.3f |\n",
            70, 30, mean(mc_results$Test_Coverage) * 100, sd(mc_results$Test_Coverage) * 100,
            mean(mc_results$Test_Width), median(mc_results$Test_Width)))
cat("\n===============================================================================\n")
cat("Execution complete. All artifacts generated successfully without modifying current package files.\n")
cat("===============================================================================\n")
