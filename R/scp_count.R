#' Spatially Localized Conformal Prediction for Count Data
#'
#' Constructs discrete, non-negative, model-agnostic prediction intervals for
#' spatially and spatio-temporally dependent count responses (e.g. disease incidences,
#' mortality cases, crime counts, or ecological abundances) using localized conformal
#' calibration with Pearson dispersion-scaled scoring.
#'
#' @param s_train Numeric matrix or data frame of training coordinates (N_train x 2).
#' @param y_train Integer or numeric vector of non-negative integer counts (length N_train).
#' @param s0 Numeric matrix or data frame of target coordinates for prediction (N_test x 2).
#' @param pred_fun Function with signature \code{function(s_train, y_train, s_new)} that
#'   fits a count model (e.g. Poisson GLM/GAM, Negative Binomial, or Random Forest) and
#'   returns a numeric vector of positive predicted mean counts for \code{s_new}.
#' @param alpha Miscoverage level; prediction intervals target \eqn{1 - \alpha}
#'   spatial coverage. Default: 0.10 (90 percent confidence).
#' @param split Proportion of training data used for model fitting; the remainder
#'   is used as the calibration set for conformal calibration. Default: 0.50.
#' @param bandwidth Optional numeric spatial kernel bandwidth. Defaults to the median
#'   pairwise distance of calibration points.
#' @param score_type Scoring method: \code{"pearson"} (standardized Pearson residual
#'   accounting for overdispersion, default), \code{"anscombe"} (variance-stabilizing
#'   square-root transform), or \code{"absolute"} (raw difference).
#'
#' @return An object of class \code{"spconform_count"} (inheriting from \code{"spconform"}) containing:
#' \describe{
#'   \item{pred}{Predicted mean counts at target locations.}
#'   \item{lower}{Integer lower prediction bounds (\eqn{\ge 0}).}
#'   \item{upper}{Integer upper prediction bounds.}
#'   \item{alpha}{Nominal miscoverage rate.}
#'   \item{bandwidth}{Spatial kernel bandwidth used.}
#'   \item{dispersion}{Estimated overdispersion parameter phi.}
#'   \item{score_type}{Scoring method used.}
#'   \item{s0}{Target spatial coordinates matrix.}
#'   \item{type}{Character string \code{"count"}.}
#' }
#' @export
#'
#' @examples
#' set.seed(123)
#' n_tr <- 100
#' s_tr <- matrix(runif(n_tr * 2), n_tr, 2)
#' lambda_tr <- exp(1.5 + 2 * s_tr[, 1] - s_tr[, 2])
#' y_tr <- rpois(n_tr, lambda = lambda_tr)
#' s_te <- matrix(runif(30 * 2), 30, 2)
#' 
#' # Define Poisson GLM predictor function
#' pred_fun <- function(s_train, y_train, s_new) {
#'   df_tr <- data.frame(y = y_train, x = s_train[, 1], y_coord = s_train[, 2])
#'   df_nw <- data.frame(x = s_new[, 1], y_coord = s_new[, 2])
#'   fit <- glm(y ~ x + y_coord, data = df_tr, family = poisson)
#'   pmax(1e-4, as.numeric(predict(fit, newdata = df_nw, type = "response")))
#' }
#' 
#' fit_cnt <- scp_count(s_tr, y_tr, s_te, pred_fun, alpha = 0.10)
#' print(fit_cnt)
scp_count <- function(s_train, y_train, s0, pred_fun, 
                      alpha = 0.10, split = 0.50, 
                      bandwidth = NULL,
                      score_type = c("pearson", "anscombe", "absolute")) {
  score_type <- match.arg(score_type)
  cl <- match.call()
  
  if (alpha <= 0 || alpha >= 1) stop("'alpha' must be strictly between 0 and 1.")
  if (split <= 0 || split >= 1) stop("'split' must be strictly between 0 and 1.")
  if (!is.function(pred_fun)) stop("'pred_fun' must be a valid prediction function.")
  
  s_train <- as.matrix(s_train)
  s0      <- as.matrix(s0)
  n_total <- nrow(s_train)
  
  if (n_total < 10) stop("Training data requires at least n = 10 observations.")
  if (length(y_train) != n_total) stop("length(y_train) must equal nrow(s_train).")
  
  y_train <- as.numeric(y_train)
  if (any(y_train < 0) || any(is.na(y_train))) {
    stop("'y_train' must contain non-negative, non-missing integer count values.")
  }
  
  # Split training data into proper training and calibration subsets
  n_fit <- max(2, round(n_total * split))
  idx_fit <- sample.int(n_total, size = n_fit)
  idx_cal <- setdiff(seq_len(n_total), idx_fit)
  n_cal <- length(idx_cal)
  
  s_fit <- s_train[idx_fit, , drop = FALSE]
  y_fit <- y_train[idx_fit]
  s_cal <- s_train[idx_cal, , drop = FALSE]
  y_cal <- y_train[idx_cal]
  
  # Fit base model and predict mean counts
  mu_cal <- pmax(1e-6, as.numeric(pred_fun(s_fit, y_fit, s_cal)))
  mu_0   <- pmax(1e-6, as.numeric(pred_fun(s_fit, y_fit, s0)))
  
  # Estimate overdispersion parameter phi: Var(Y) = mu + phi * mu^2
  # Method of moments on calibration residuals
  sq_res <- (y_cal - mu_cal)^2
  phi_hat <- max(0, mean((sq_res - mu_cal) / (mu_cal^2), na.rm = TRUE))
  if (!is.finite(phi_hat)) phi_hat <- 0
  
  # Evaluate nonconformity scores on calibration set
  if (score_type == "pearson") {
    # Standardized Pearson residual: |y - mu| / sqrt(mu + phi * mu^2)
    scale_cal <- sqrt(mu_cal + phi_hat * (mu_cal^2))
    cal_scores <- abs(y_cal - mu_cal) / scale_cal
  } else if (score_type == "anscombe") {
    # Anscombe variance-stabilizing transform: 2 * |sqrt(y + 3/8) - sqrt(mu + 3/8)|
    cal_scores <- 2 * abs(sqrt(y_cal + 3/8) - sqrt(mu_cal + 3/8))
  } else {
    # Raw absolute differences
    cal_scores <- abs(y_cal - mu_cal)
  }
  
  # Determine spatial bandwidth h
  if (is.null(bandwidth) || bandwidth <= 0) {
    dist_cal <- as.matrix(stats::dist(s_cal))
    d_vals <- dist_cal[lower.tri(dist_cal)]
    bandwidth <- stats::median(d_vals[d_vals > 0])
  }
  
  # Distance matrix from calibration points to target locations: (n_cal x n_test)
  diff_x <- outer(s_cal[, 1], s0[, 1], "-")
  diff_y <- outer(s_cal[, 2], s0[, 2], "-")
  dist_cal_0 <- sqrt(diff_x^2 + diff_y^2)
  
  n_test <- nrow(s0)
  lower <- integer(n_test)
  upper <- integer(n_test)
  
  for (j in seq_len(n_test)) {
    # Spatial kernel weights: w_i(s0_j)
    w_vec <- exp(- (dist_cal_0[, j]^2) / (2 * bandwidth^2)) + 1e-12
    w_norm <- w_vec / sum(w_vec)
    
    ord_sc <- order(cal_scores)
    scores_sorted <- cal_scores[ord_sc]
    w_sorted <- w_norm[ord_sc]
    cum_w <- cumsum(w_sorted)
    
    q_level <- min(1.0, (1 - alpha) * (1 + 1 / n_cal))
    idx_q <- which(cum_w >= q_level)[1]
    if (is.na(idx_q)) idx_q <- length(scores_sorted)
    q_val <- scores_sorted[idx_q]
    
    mu_j <- mu_0[j]
    
    if (score_type == "pearson") {
      scale_j <- sqrt(mu_j + phi_hat * (mu_j^2))
      raw_lwr <- mu_j - q_val * scale_j
      raw_upr <- mu_j + q_val * scale_j
    } else if (score_type == "anscombe") {
      t_mu <- sqrt(mu_j + 3/8)
      raw_lwr <- max(0, (t_mu - q_val / 2)^2 - 3/8)
      raw_upr <- (t_mu + q_val / 2)^2 - 3/8
    } else {
      raw_lwr <- mu_j - q_val
      raw_upr <- mu_j + q_val
    }
    
    # Project onto strictly non-negative integer domain
    lower[j] <- max(0L, as.integer(floor(raw_lwr)))
    upper[j] <- max(lower[j], as.integer(ceiling(raw_upr)))
  }
  
  structure(
    list(
      pred = mu_0,
      lower = lower,
      upper = upper,
      alpha = alpha,
      s0 = s0,
      bandwidth = round(bandwidth, 4),
      dispersion = round(phi_hat, 4),
      score_type = score_type,
      type = "count",
      call = cl
    ),
    class = c("spconform_count", "spconform")
  )
}

#' @export
print.spconform_count <- function(x, ...) {
  cat("-----------------------------------------------------------\n")
  cat("  <spconform> Spatial Conformal Count Prediction\n")
  cat("-----------------------------------------------------------\n")
  cat(sprintf("Target Nominal Coverage:      %.1f%%\n", 100 * (1 - x$alpha)))
  cat("Score Formulation:            ", x$score_type, "\n", sep = "")
  cat(sprintf("Number of Prediction Points:  %d\n", length(x$pred)))
  cat(sprintf("Spatial Bandwidth (h):        %.4f\n", x$bandwidth))
  cat(sprintf("Overdispersion (phi):         %.4f (%s)\n", 
              x$dispersion, if (x$dispersion > 0.05) "Overdispersed / NegBin" else "Equidispersed / Poisson"))
  cat("-----------------------------------------------------------\n")
  widths <- x$upper - x$lower
  cat(sprintf("Predicted Count Range:        [%.2f, %.2f]\n", min(x$pred), max(x$pred)))
  cat(sprintf("Mean Integer Interval Width:  %.2f counts\n", mean(widths)))
  cat(sprintf("Zero Lower-Bound Locations:   %d (%.1f%%)\n", sum(x$lower == 0), 100 * mean(x$lower == 0)))
  cat("-----------------------------------------------------------\n")
  invisible(x)
}

#' @export
summary.spconform_count <- function(object, ...) {
  print(object)
  widths <- object$upper - object$lower
  cat("\nInteger Prediction Interval Width Quantiles:\n")
  print(stats::quantile(widths, probs = c(0, 0.25, 0.5, 0.75, 1.0)))
  invisible(object)
}

#' @export
plot.spconform_count <- function(x, ...) {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  
  graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  
  widths <- x$upper - x$lower
  # Color scale based on predicted count
  val_scaled <- (x$pred - min(x$pred)) / (max(x$pred) - min(x$pred) + 1e-12)
  pal <- grDevices::colorRampPalette(c("#4575B4", "#FFFFBF", "#D73027"))(100)
  cols <- pal[pmax(1, pmin(100, round(val_scaled * 99) + 1))]
  
  # Panel 1: Geographic Map of Predicted Counts with Width Bubble Size
  cex_size <- 0.8 + 1.2 * (widths / max(widths))
  plot(x$s0[, 1], x$s0[, 2], col = cols, pch = 19, cex = cex_size,
       xlab = "Easting (X)", ylab = "Northing (Y)",
       main = "Predicted Counts & Interval Widths")
  graphics::grid(col = "gray85")
  graphics::points(x$s0[, 1], x$s0[, 2], col = "gray30", pch = 1, cex = cex_size)
  
  # Panel 2: Predicted Counts vs Integer Prediction Bounds
  ord <- order(x$pred)
  idx_seq <- seq_along(ord)
  plot(idx_seq, x$pred[ord], type = "n",
       ylim = c(0, max(x$upper[ord]) * 1.1),
       xlab = "Locations (Ranked by Predicted Count)",
       ylab = "Count (Integers >= 0)",
       main = "Discrete Prediction Intervals")
  graphics::grid(col = "gray85")
  
  # Shaded interval area
  graphics::polygon(c(idx_seq, rev(idx_seq)),
                    c(x$lower[ord], rev(x$upper[ord])),
                    col = grDevices::adjustcolor("#4575B4", alpha.f = 0.3),
                    border = NA)
  graphics::lines(idx_seq, x$lower[ord], col = "#4575B4", lty = 2, lwd = 1.2)
  graphics::lines(idx_seq, x$upper[ord], col = "#4575B4", lty = 2, lwd = 1.2)
  graphics::lines(idx_seq, x$pred[ord], col = "#D73027", lwd = 2)
  graphics::abline(h = 0, lty = 1, col = "black")
  
  graphics::legend("topleft", legend = c("Predicted Mean", "Conformal Bounds (>= 0)"),
                   col = c("#D73027", "#4575B4"), lwd = c(2, 1.5), lty = c(1, 2),
                   bty = "n", cex = 0.8)
}
