#' Spatial Conditional Coverage Test for Conformal Prediction
#'
#' Formally tests whether the coverage of conformal prediction intervals varies
#' systematically across spatial coordinates using a spatial logistic likelihood
#' ratio test. Evaluates the null hypothesis H0: P(Y(s) in C(s) | S = s) = 1 - alpha.
#' Includes automatic degree regularization and finite-sample permutation testing
#' to prevent quasi-complete separation artifacts when miscoverage counts are small.
#'
#' @param object An object of class \code{"spconform"}.
#' @param y_true Numeric vector of observed test responses.
#' @param coords Optional numeric matrix or data frame of spatial coordinates (N x 2).
#'   If \code{NULL}, coordinates are extracted from \code{object$coords}.
#' @param degree Integer degree of the spatial polynomial trend (default: 2, quadratic).
#'   Automatically reduced to 1 if the number of miscovered points is small to prevent
#'   overfitting and separation.
#' @param test Character specifying the test statistic: \code{"LRT"} (Likelihood
#'   Ratio Test, default) or \code{"Wald"}.
#' @param permutation Logical or \code{"auto"}. If \code{"auto"} (default), runs a
#'   finite-sample spatial permutation test whenever miscoverage count <= 10 or N <= 150.
#' @param n_perm Number of permutation replicates when permutation testing is active (default: 499).
#'
#' @return An object of class \code{"spconform_coverage_test"} containing:
#' \describe{
#'   \item{statistic}{The chi-squared test statistic value.}
#'   \item{df}{Degrees of freedom under the null distribution.}
#'   \item{p_value}{Reported p-value (permutation p-value if permutation active, otherwise asymptotic).}
#'   \item{p_asymptotic}{Asymptotic chi-squared p-value.}
#'   \item{p_perm}{Permutation test p-value (or NA if not computed).}
#'   \item{empirical_coverage}{Observed marginal coverage rate.}
#'   \item{nominal_coverage}{Target nominal coverage 1 - alpha.}
#'   \item{decision}{Verbal interpretation of the hypothesis test.}
#'   \item{hit_indicator}{Binary integer vector of coverage indicators (1 = covered, 0 = missed).}
#'   \item{method}{Name of the test performed.}
#'   \item{degree_used}{Effective polynomial degree evaluated.}
#' }
#' @export
#'
#' @examples
#' set.seed(123)
#' s_train <- matrix(runif(100 * 2), 100, 2)
#' y_train <- rnorm(100)
#' s_test  <- matrix(runif(40 * 2), 40, 2)
#' y_test  <- rnorm(40)
#' fit <- scp_geostatistical(s_train, y_train, s_test, function(tr_s, tr_y, nw_s) rep(0, nrow(nw_s)), alpha = 0.1)
#' test_res <- spatial_conditional_coverage(fit, y_test, s_test)
#' print(test_res)
spatial_conditional_coverage <- function(object, y_true, coords = NULL, degree = 2, 
                                         test = c("LRT", "Wald"), 
                                         permutation = "auto", n_perm = 499) {
  test <- match.arg(test)
  
  if (!inherits(object, "spconform")) {
    stop("'object' must be an object of class 'spconform'.")
  }
  
  n <- length(object$pred)
  if (length(y_true) != n) {
    stop(sprintf("length(y_true) (%d) must match length(object$pred) (%d).", length(y_true), n))
  }
  
  if (is.null(coords)) {
    if (!is.null(object$s0)) {
      coords <- object$s0
    } else if (!is.null(object$coords)) {
      coords <- object$coords
    } else {
      stop("Spatial coordinates must be provided via 'coords' or present in 'object$s0'.")
    }
  }
  
  coords <- as.matrix(coords)
  if (nrow(coords) != n || ncol(coords) < 2) {
    stop("'coords' must have at least 2 columns and match the number of predictions.")
  }
  
  # Binary coverage hit indicator: 1 = covered, 0 = miscovered
  hit <- as.integer(y_true >= object$lower & y_true <= object$upper)
  n_hit  <- sum(hit == 1)
  n_miss <- sum(hit == 0)
  emp_cov <- mean(hit)
  nom_cov <- 1 - object$alpha
  
  # Edge case: perfect homogeneous coverage (all 1s or all 0s)
  if (n_miss == 0 || n_hit == 0) {
    res <- list(
      statistic = 0,
      df = 0,
      p_value = 1.0,
      p_asymptotic = 1.0,
      p_perm = 1.0,
      empirical_coverage = emp_cov,
      nominal_coverage = nom_cov,
      decision = "Uniform coverage achieved (no spatial variance to test).",
      hit_indicator = hit,
      method = paste("Spatial Conditional Coverage (Homogeneous)", test),
      degree_used = 0,
      model = NULL
    )
    class(res) <- "spconform_coverage_test"
    return(res)
  }
  
  # Guard against quasi-complete separation:
  # If minority class has <= 4 observations, reduce degree to 1 (x + y only, df = 2)
  min_class <- min(n_hit, n_miss)
  degree_used <- degree
  if (min_class <= 4 && degree > 1) {
    degree_used <- 1
  }
  
  # Center coordinates for numerical stability
  x_coord <- scale(coords[, 1], center = TRUE, scale = TRUE)
  y_coord <- scale(coords[, 2], center = TRUE, scale = TRUE)
  
  df_model <- data.frame(hit = hit, x = as.numeric(x_coord), y = as.numeric(y_coord))
  
  # Null model: spatially invariant coverage logit(P(Hit=1)) = beta_0
  m_null <- suppressWarnings(
    stats::glm(hit ~ 1, data = df_model, family = stats::binomial(link = "logit"))
  )
  
  # Spatial alternative model: coverage depends on spatial coordinates
  if (degree_used == 1) {
    m_spatial <- suppressWarnings(
      stats::glm(hit ~ x + y, data = df_model, family = stats::binomial(link = "logit"))
    )
  } else {
    m_spatial <- suppressWarnings(
      stats::glm(hit ~ stats::poly(x, y, degree = degree_used, raw = TRUE), 
                 data = df_model, family = stats::binomial(link = "logit"))
    )
  }
  
  if (test == "LRT") {
    # Likelihood Ratio Test: deviance difference
    stat <- stats::deviance(m_null) - stats::deviance(m_spatial)
    stat <- max(0, stat)
    df_diff <- m_null$df.residual - m_spatial$df.residual
    p_asymp <- stats::pchisq(stat, df = df_diff, lower.tail = FALSE)
  } else {
    # Wald Test on all spatial coefficients
    coefs <- stats::coef(m_spatial)[-1] # drop intercept
    vcov_mat <- stats::vcov(m_spatial)[-1, -1, drop = FALSE]
        base_pinv <- function(A, tol = 1e-7) {
      s <- svd(A)
      pos <- s$d > tol * max(s$d)
      if (all(!pos)) return(matrix(0, ncol(A), nrow(A)))
      s$v[, pos, drop = FALSE] %*% ((1 / s$d[pos]) * t(s$u[, pos, drop = FALSE]))
    }
    inv_vcov <- tryCatch(solve(vcov_mat), error = function(e) base_pinv(vcov_mat))
    stat <- as.numeric(t(coefs) %*% inv_vcov %*% coefs)
    stat <- max(0, stat)
    df_diff <- length(coefs)
    p_asymp <- stats::pchisq(stat, df = df_diff, lower.tail = FALSE)
  }
  
  # Decide whether to run permutation test
  do_perm <- if (identical(permutation, "auto")) {
    (min_class <= 10 || n <= 150)
  } else {
    isTRUE(permutation)
  }
  
  p_perm <- NA_real_
  if (do_perm && n_perm > 0) {
    perm_stats <- numeric(n_perm)
    for (b in seq_len(n_perm)) {
      df_perm <- df_model
      df_perm$hit <- sample(df_model$hit)
      
      m_null_b <- suppressWarnings(
        stats::glm(hit ~ 1, data = df_perm, family = stats::binomial(link = "logit"))
      )
      if (degree_used == 1) {
        m_sp_b <- suppressWarnings(
          stats::glm(hit ~ x + y, data = df_perm, family = stats::binomial(link = "logit"))
        )
      } else {
        m_sp_b <- suppressWarnings(
          stats::glm(hit ~ stats::poly(x, y, degree = degree_used, raw = TRUE), 
                     data = df_perm, family = stats::binomial(link = "logit"))
        )
      }
      
      dev_diff <- stats::deviance(m_null_b) - stats::deviance(m_sp_b)
      perm_stats[b] <- max(0, dev_diff)
    }
    # Finite-sample exact permutation p-value
    p_perm <- (1 + sum(perm_stats >= stat)) / (n_perm + 1)
  }
  
  # Active reported p-value
  final_p <- if (!is.na(p_perm)) p_perm else p_asymp
  
  decision <- if (final_p >= 0.05) {
    "Spatial conditional coverage holds (Fail to reject H0; no geographic bias)."
  } else {
    "Spatial conditional coverage violated (Reject H0; coverage varies across space)."
  }
  
  res <- list(
    statistic = round(stat, 4),
    df = df_diff,
    p_value = final_p,
    p_asymptotic = p_asymp,
    p_perm = p_perm,
    empirical_coverage = emp_cov,
    nominal_coverage = nom_cov,
    decision = decision,
    hit_indicator = hit,
    method = paste0("Spatial Conditional Coverage (", test, ", degree = ", degree_used, 
                    if (!is.na(p_perm)) paste0(", B = ", n_perm, " perms)") else ")"),
    degree_used = degree_used,
    model = m_spatial
  )
  class(res) <- "spconform_coverage_test"
  return(res)
}

#' @export
print.spconform_coverage_test <- function(x, ...) {
  cat("-----------------------------------------------------------\n")
  cat("  spconform: Spatial Conditional Coverage Test\n")
  cat("-----------------------------------------------------------\n")
  cat("Test Method:           ", x$method, "\n", sep = "")
  cat(sprintf("Nominal Target:        %.1f%%\n", x$nominal_coverage * 100))
  cat(sprintf("Empirical Coverage:    %.1f%%\n", x$empirical_coverage * 100))
  cat(sprintf("Test Statistic (Chi2): %.4f (df = %d)\n", x$statistic, x$df))
  cat("p-value:               ", format.pval(x$p_value, eps = 1e-4, digits = 4), 
      if (!is.na(x$p_perm)) " (permutation exact)" else " (asymptotic)", "\n", sep = "")
  if (!is.na(x$p_perm)) {
    cat("  - Asymptotic p-value: ", format.pval(x$p_asymptotic, eps = 1e-4, digits = 4), "\n", sep = "")
    cat("  - Permutation p-value: ", format.pval(x$p_perm, eps = 1e-4, digits = 4), "\n", sep = "")
  }
  cat("Status:                ", if (x$p_value >= 0.05) "PASS (H0 accepted)" else "FAIL (H0 rejected)", "\n", sep = "")
  cat("Conclusion:            ", x$decision, "\n", sep = "")
  cat("-----------------------------------------------------------\n")
  invisible(x)
}


#' Geary's C Test for Spatial Autocorrelation
#'
#' Computes Geary's C statistic to test for local spatial autocorrelation and
#' spatial variance discontinuities in conformal prediction residuals or interval widths.
#' Follows Cliff and Ord (1981) analytical moments under randomization.
#'
#' @param x Numeric vector of values (e.g. residuals, interval widths, or nonconformity scores).
#' @param coords Numeric matrix or data frame of spatial coordinates (N x 2).
#' @param bandwidth Optional numeric spatial bandwidth. Defaults to the median pairwise distance.
#' @param alternative Character string specifying the alternative hypothesis:
#'   \code{"two.sided"} (default), \code{"greater"} (positive spatial autocorrelation, C < 1),
#'   or \code{"less"} (negative spatial autocorrelation / spatial discontinuity, C > 1).
#'
#' @return An object of class \code{"spconform_geary"} containing:
#' \describe{
#'   \item{C}{Observed Geary's C statistic.}
#'   \item{expected}{Expected value under the randomization null hypothesis (1.0).}
#'   \item{variance}{Analytical variance of C under randomization.}
#'   \item{z_score}{Standardized test statistic.}
#'   \item{p_value}{Asymptotic p-value.}
#'   \item{alternative}{The alternative hypothesis specified.}
#' }
#' @export
#'
#' @examples
#' set.seed(123)
#' coords <- matrix(runif(80 * 2), 80, 2)
#' resids <- rnorm(80)
#' geary_res <- spatial_geary(resids, coords)
#' print(geary_res)
spatial_geary <- function(x, coords, bandwidth = NULL, 
                          alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)
  
  x <- as.numeric(x)
  n <- length(x)
  if (n < 4) {
    stop("Geary's C test requires at least n = 4 observations.")
  }
  
  coords <- as.matrix(coords)
  if (nrow(coords) != n || ncol(coords) < 2) {
    stop("nrow(coords) must equal length(x) with at least 2 columns.")
  }
  
  # Remove missing values
  valid <- !is.na(x) & !is.na(coords[, 1]) & !is.na(coords[, 2])
  if (sum(valid) < n) {
    x <- x[valid]
    coords <- coords[valid, , drop = FALSE]
    n <- length(x)
  }
  
  # Compute pairwise Euclidean distance matrix
  dist_mat <- as.matrix(stats::dist(coords))
  
  # Spatial bandwidth h
  if (is.null(bandwidth) || bandwidth <= 0) {
    d_vals <- dist_mat[lower.tri(dist_mat)]
    bandwidth <- stats::median(d_vals[d_vals > 0])
  }
  
  # Gaussian spatial weight matrix: w_ij = exp(-d_ij^2 / (2 * h^2))
  W_mat <- exp(- (dist_mat^2) / (2 * bandwidth^2))
  diag(W_mat) <- 0
  
  S0 <- sum(W_mat)
  if (S0 == 0) {
    stop("Sum of spatial weights is zero; check coordinate scales or bandwidth.")
  }
  
  # Deviations from mean
  z <- scale(x, scale = TRUE)
  zz <- sum(z^2)
  if (zz == 0) {
    stop("Variance of x is zero; cannot compute Geary's C.")
  }
  
  # Numerator of Geary's C: sum_i sum_j w_ij (x_i - x_j)^2
  diff_mat <- outer(as.numeric(z), as.numeric(z), "-")^2
  C_obs <- ((n - 1) / (2 * S0)) * (sum(W_mat * diff_mat) / zz)
  
  # Expected value under H0
  E_C <- 1.0
  
  # Exact Analytical Variance under Randomization (Cliff and Ord 1973, 1981)
  S02 <- S0 * S0
  S1 <- 2 * sum(W_mat^2)
  S2 <- 4 * sum(rowSums(W_mat)^2)
  
  K <- (n * sum(z^4)) / (zz^2)
  n1 <- n - 1
  n2 <- n - 2
  n3 <- n - 3
  nn <- n * n
  
  VC <- (n1 * S1 * (nn - 3 * n + 3 - K * n1))
  VC <- VC - ((1 / 4) * (n1 * S2 * (nn + 3 * n - 6 - K * (nn - n + 2))))
  VC <- VC + (S02 * (nn - 3 - K * (n1^2)))
  VC <- VC / (n * n2 * n3 * S02)
  VC <- max(VC, 1e-12) # Numerical floor safeguard
  
  # Standardized test statistic:
  z_score <- (E_C - C_obs) / sqrt(VC)
  
  p_val <- switch(alternative,
    "two.sided" = 2 * stats::pnorm(abs(z_score), lower.tail = FALSE),
    "greater"   = stats::pnorm(z_score, lower.tail = FALSE), # positive spatial autocorr (C < 1)
    "less"      = stats::pnorm(z_score, lower.tail = TRUE)    # negative spatial autocorr (C > 1)
  )
  
  res <- list(
    C = round(C_obs, 4),
    expected = E_C,
    variance = round(VC, 6),
    z_score = round(z_score, 4),
    p_value = p_val,
    bandwidth = round(bandwidth, 4),
    alternative = alternative,
    sample_size = n
  )
  class(res) <- "spconform_geary"
  return(res)
}

#' @export
print.spconform_geary <- function(x, ...) {
  cat("-----------------------------------------------------------\n")
  cat("  spconform: Geary's C Spatial Autocorrelation Test\n")
  cat("-----------------------------------------------------------\n")
  cat(sprintf("Observed Geary's C:    %.4f\n", x$C))
  cat(sprintf("Expectation (H0):      %.4f\n", x$expected))
  cat(sprintf("Analytical Variance:   %.6f\n", x$variance))
  cat(sprintf("Standardized z-score:  %.4f\n", x$z_score))
  cat("Asymptotic p-value:    ", format.pval(x$p_value, eps = 1e-4, digits = 4), "\n", sep = "")
  cat(sprintf("Spatial Bandwidth:     %.4f\n", x$bandwidth))
  cat("Alternative:           ", x$alternative, " (", 
      if (x$C < 1) "positive spatial autocorrelation / clustering" else "negative spatial autocorrelation / local disparity",
      ")\n", sep = "")
  cat("Status:                ", if (x$p_value < 0.05) "SIGNIFICANT (Spatial structure detected)" else "NOT SIGNIFICANT (Spatial randomness)", "\n", sep = "")
  cat("-----------------------------------------------------------\n")
  invisible(x)
}


#' Local Indicators of Spatial Association (Anselin's LISA / Local Moran's I)
#'
#' Computes Local Moran's I_i statistics for each observation to detect spatial clusters
#' (hotspots/coldspots) and spatial outliers of conformal prediction uncertainty,
#' following Luc Anselin (1995). Supports direct diagnostics of prediction residuals,
#' interval widths (uncertainty zoning), or raw numeric vectors, with multiple testing
#' adjustments (FDR / Bonferroni).
#'
#' @param object An object of class \code{"spconform"} or a numeric vector.
#' @param y_true Optional numeric vector of true test responses (required when
#'   \code{object} is of class \code{"spconform"} and \code{target = "residuals"}).
#' @param coords Optional numeric matrix or data frame of spatial coordinates (N x 2).
#'   Extracted automatically if present in \code{object$coords}.
#' @param target Character string specifying the diagnostic target when \code{object}
#'   is of class \code{"spconform"}: \code{"residuals"} (absolute errors |y - y_hat| for
#'   identifying systematic local model failure) or \code{"widths"} (interval widths for
#'   spatial uncertainty zoning).
#' @param bandwidth Optional numeric spatial bandwidth. Defaults to the median pairwise distance.
#' @param p_adjust_method Method for p-value adjustment across the N multiple hypothesis
#'   tests (passed to \code{stats::p.adjust}): \code{"fdr"} (default, Benjamini & Hochberg),
#'   \code{"bonferroni"}, \code{"holm"}, or \code{"none"}.
#' @param significance_level Alpha threshold for declaring a local cluster significant
#'   (default: 0.05).
#'
#' @return An object of class \code{"spconform_lisa"} containing:
#' \describe{
#'   \item{Ii}{Numeric vector of Local Moran's I_i statistics.}
#'   \item{expected}{Expected value under conditional randomization.}
#'   \item{variance}{Analytical variance of I_i under randomization.}
#'   \item{z_score}{Standardized z-scores for each location.}
#'   \item{p_raw}{Raw asymptotic two-sided p-values.}
#'   \item{p_adjusted}{Adjusted p-values under the chosen multiple testing method.}
#'   \item{quadrant}{Character vector classifying each location into:
#'     \code{"High-High"}, \code{"Low-Low"}, \code{"High-Low"}, \code{"Low-High"},
#'     or \code{"Not Significant"}.}
#'   \item{spatial_lag}{Spatially lagged values of centered x.}
#'   \item{coords}{Spatial coordinates matrix.}
#'   \item{bandwidth}{Spatial bandwidth used.}
#'   \item{target}{The diagnostic target evaluated.}
#'   \item{p_adjust_method}{Multiple testing correction method applied.}
#' }
#' @export
#'
#' @examples
#' set.seed(123)
#' coords <- matrix(runif(60 * 2), 60, 2)
#' resids <- rnorm(60)
#' lisa_res <- spatial_lisa(resids, coords = coords)
#' print(lisa_res)
spatial_lisa <- function(object, coords = NULL, y_true = NULL,
                         target = c("residuals", "widths"),
                         bandwidth = NULL, 
                         p_adjust_method = c("fdr", "bonferroni", "holm", "none"),
                         significance_level = 0.05) {
  target <- match.arg(target)
  p_adjust_method <- match.arg(p_adjust_method)
  
  if (inherits(object, "spconform")) {
    if (is.null(coords)) {
      if (!is.null(object$s0)) {
        coords <- object$s0
      } else if (!is.null(object$coords)) {
        coords <- object$coords
      } else {
        stop("Spatial coordinates must be provided via 'coords' or present in 'object$s0'.")
      }
    }
    
    if (target == "residuals") {
      if (is.null(y_true)) {
        stop("'y_true' must be provided when target = 'residuals'.")
      }
      if (length(y_true) != length(object$pred)) {
        stop("length(y_true) must match length(object$pred).")
      }
      x <- abs(y_true - object$pred) # Absolute prediction errors
      target_label <- "Absolute Residuals |Y - Y_hat| (Model Bias Diagnostic)"
    } else {
      x <- object$upper - object$lower # Prediction interval widths
      target_label <- "Interval Widths (Spatial Uncertainty Zoning)"
    }
  } else if (is.numeric(object)) {
    x <- as.numeric(object)
    if (is.null(coords)) {
      stop("When 'object' is a numeric vector, 'coords' must be provided.")
    }
    target_label <- "Numeric Vector"
  } else {
    stop("'object' must be an object of class 'spconform' or a numeric vector.")
  }
  
  n <- length(x)
  if (n < 5) {
    stop("LISA test requires at least n = 5 observations.")
  }
  
  coords <- as.matrix(coords)
  if (nrow(coords) != n || ncol(coords) < 2) {
    stop("nrow(coords) must equal length of data with at least 2 coordinate columns.")
  }
  
  # Remove missing values
  valid <- !is.na(x) & !is.na(coords[, 1]) & !is.na(coords[, 2])
  if (sum(valid) < n) {
    x <- x[valid]
    coords <- coords[valid, , drop = FALSE]
    n <- length(x)
  }
  
  # Compute pairwise Euclidean distance matrix
  dist_mat <- as.matrix(stats::dist(coords))
  
  # Spatial bandwidth h
  if (is.null(bandwidth) || bandwidth <= 0) {
    d_vals <- dist_mat[lower.tri(dist_mat)]
    bandwidth <- stats::median(d_vals[d_vals > 0])
  }
  
  # Gaussian spatial weight matrix: w_ij = exp(-d_ij^2 / (2 * h^2))
  W_mat <- exp(- (dist_mat^2) / (2 * bandwidth^2))
  diag(W_mat) <- 0 # Strict zero diagonal
  
  # Row-standardize weights so sum_j w_ij = 1
  row_sums <- rowSums(W_mat)
  row_sums[row_sums == 0] <- 1
  W_std <- W_mat / row_sums
  
  # Centered variable z_i = x_i - mean(x)
  x_bar <- mean(x)
  z <- x - x_bar
  m2 <- sum(z^2) / n
  if (m2 == 0) {
    stop("Variance of diagnostic variable is zero; cannot compute LISA.")
  }
  
  # Spatial lag: lag_i = sum_j w_std_ij * z_j
  spatial_lag <- as.numeric(W_std %*% z)
  
  # Local Moran's I_i = (z_i / m2) * lag_i
  I_i <- (z / m2) * spatial_lag
  
  # Sample kurtosis b2 = (1/n * sum(z^4)) / (m2^2)
  b2 <- (sum(z^4) / n) / (m2^2)
  
  # Analytical Conditional Variance (Anselin 1995, eq. 8-10):
  E_Ii <- - 1 / (n - 1)
  w_i2 <- rowSums(W_std^2)
  
  term1 <- w_i2 * (n - b2) / (n - 1)
  term2 <- (1 - w_i2) * (2 * b2 - n) / ((n - 1) * (n - 2))
  var_Ii <- term1 + term2 - (E_Ii^2)
  var_Ii <- pmax(var_Ii, 1e-12) # Floor safeguard
  
  # Standardized z-score
  z_score <- (I_i - E_Ii) / sqrt(var_Ii)
  
  # Asymptotic two-sided p-values
  p_raw <- 2 * stats::pnorm(abs(z_score), lower.tail = FALSE)
  
  # Multiple testing correction
  p_adj <- stats::p.adjust(p_raw, method = p_adjust_method)
  
  # Quadrant classification:
  quad_raw <- ifelse(z >= 0 & spatial_lag >= 0, "High-High",
              ifelse(z < 0 & spatial_lag < 0, "Low-Low",
              ifelse(z >= 0 & spatial_lag < 0, "High-Low", "Low-High")))
  
  quadrant <- ifelse(p_adj < significance_level, quad_raw, "Not Significant")
  quadrant <- factor(quadrant, levels = c("High-High", "Low-Low", "High-Low", "Low-High", "Not Significant"))
  
  res <- list(
    Ii = round(I_i, 4),
    expected = round(E_Ii, 4),
    variance = round(var_Ii, 6),
    z_score = round(z_score, 4),
    p_raw = p_raw,
    p_adjusted = p_adj,
    quadrant = quadrant,
    spatial_lag = round(spatial_lag, 4),
    z = round(z, 4),
    coords = coords,
    bandwidth = round(bandwidth, 4),
    target = target_label,
    p_adjust_method = p_adjust_method,
    significance_level = significance_level,
    n = n
  )
  class(res) <- "spconform_lisa"
  return(res)
}

#' @export
print.spconform_lisa <- function(x, ...) {
  cat("-----------------------------------------------------------\n")
  cat("  spconform: Anselin's Local Moran's I (LISA) Test\n")
  cat("-----------------------------------------------------------\n")
  cat("Diagnostic Target:            ", x$target, "\n", sep = "")
  cat(sprintf("Number of Spatial Locations:  %d\n", x$n))
  cat(sprintf("Spatial Bandwidth (h):        %.4f\n", x$bandwidth))
  cat("Multiple Testing Adjustment:  ", x$p_adjust_method, "\n", sep = "")
  cat(sprintf("Significance Alpha Level:     %.2f\n", x$significance_level))
  cat("-----------------------------------------------------------\n")
  cat("LISA Cluster Classification Breakdown:\n")
  counts <- table(x$quadrant)
  pcts <- prop.table(counts) * 100
  for (nm in names(counts)) {
    cat(sprintf("  %-16s : %3d (%5.1f%%)\n", nm, counts[nm], pcts[nm]))
  }
  cat("-----------------------------------------------------------\n")
  invisible(x)
}

#' @export
plot.spconform_lisa <- function(x, ...) {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  
  graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  
  q_colors <- c(
    "High-High"        = "#D73027", # Deep Red (Hotspot)
    "Low-Low"          = "#4575B4", # Deep Blue (Coldspot)
    "High-Low"         = "#FDAE61", # Orange (Spatial Outlier)
    "Low-High"         = "#ABD9E9", # Light Blue (Spatial Outlier)
    "Not Significant"  = "#E0E0E0"  # Neutral Light Gray
  )
  
  pt_cols <- q_colors[as.character(x$quadrant)]
  
  # Panel 1: LISA Cluster Map
  plot(x$coords[, 1], x$coords[, 2], col = pt_cols, pch = 19, cex = 1.4,
       xlab = "Easting (X)", ylab = "Northing (Y)",
       main = paste("LISA Map:", if (grepl("Residuals", x$target)) "Model Residuals" else "Uncertainty Zoning"))
  graphics::grid(col = "gray85")
  graphics::points(x$coords[, 1], x$coords[, 2], col = "gray30", pch = 1, cex = 1.4)
  
  graphics::legend("topleft", legend = names(q_colors), col = q_colors,
                   pch = 19, pt.cex = 1.2, bty = "n", cex = 0.8)
  
  # Panel 2: Moran Scatterplot (z vs Spatial Lag)
  plot(x$z, x$spatial_lag, col = pt_cols, pch = 19, cex = 1.3,
       xlab = "Standardized Value (z)", ylab = "Spatial Lag Wz",
       main = "Moran Scatterplot")
  graphics::grid(col = "gray85")
  graphics::abline(h = 0, v = 0, lty = 2, col = "gray50")
  fit_moran <- stats::lm(x$spatial_lag ~ x$z)
  graphics::abline(fit_moran, col = "black", lwd = 2)
  graphics::points(x$z, x$spatial_lag, col = "gray30", pch = 1, cex = 1.3)
}

#' Directional Anisotropy Test for Spatial Dependence
#'
#' Tests whether spatial autocorrelation in conformal prediction residuals or interval
#' widths exhibits directional anisotropy (i.e. spatial dependence varying with geographic
#' orientation, such as river corridors or prevailing winds). Formally tests H0 (Isotropy)
#' versus H1 (Anisotropy) using directional semivariance ANOVA and calculates the
#' principal axis of spatial continuity and the geometric anisotropy ratio.
#'
#' @param object An object of class \code{"spconform"} or a numeric vector of values.
#' @param coords Optional numeric matrix or data frame of spatial coordinates (N x 2).
#'   Extracted automatically if present in \code{object$s0} or \code{object$coords}.
#' @param y_true Optional numeric vector of true responses (used when \code{target = "residuals"}).
#' @param target Diagnostic target: \code{"residuals"} (default, |y - y_hat|) or \code{"widths"}.
#' @param n_directions Integer number of directional bins (default: 4, corresponding to
#'   0 deg [E-W], 45 deg [NE-SW], 90 deg [N-S], and 135 deg [NW-SE]).
#' @param max_dist Maximum pairwise distance lag considered (default: 0.5 * max pairwise distance).
#' @param n_perm Number of spatial coordinate rotation permutations for finite-sample testing (default: 299).
#'
#' @return An object of class \code{"spconform_anisotropy"} containing:
#' \describe{
#'   \item{F_statistic}{The ANOVA F test statistic testing directional heterogeneity.}
#'   \item{p_value}{Reported p-value (permutation exact if active, otherwise asymptotic).}
#'   \item{p_asymptotic}{Asymptotic F-distribution p-value.}
#'   \item{p_perm}{Permutation p-value.}
#'   \item{anisotropy_ratio}{Ratio of maximum to minimum directional semivariance.}
#'   \item{principal_angle}{Estimated angle (degrees) of maximum spatial continuity.}
#'   \item{directional_table}{Summary data frame of directional semivariances and sample pair counts.}
#'   \item{decision}{Verbal interpretation of the isotropy/anisotropy hypothesis.}
#'   \item{target}{Target variable analyzed.}
#' }
#' @export
#'
#' @examples
#' set.seed(123)
#' coords <- matrix(runif(60 * 2), 60, 2)
#' resids <- rnorm(60)
#' aniso_res <- spatial_anisotropy(resids, coords = coords)
#' print(aniso_res)
spatial_anisotropy <- function(object, coords = NULL, y_true = NULL,
                               target = c("residuals", "widths"),
                               n_directions = 4,
                               max_dist = NULL,
                               n_perm = 299) {
  target <- match.arg(target)
  
  if (inherits(object, "spconform")) {
    if (is.null(coords)) {
      if (!is.null(object$s0)) {
        coords <- object$s0
      } else if (!is.null(object$coords)) {
        coords <- object$coords
      } else {
        stop("Spatial coordinates must be provided via 'coords' or present in 'object$s0'.")
      }
    }
    if (target == "residuals") {
      if (is.null(y_true)) {
        stop("'y_true' must be provided when target = 'residuals'.")
      }
      x <- abs(y_true - object$pred)
      target_label <- "Absolute Residuals |Y - Y_hat|"
    } else {
      x <- object$upper - object$lower
      target_label <- "Prediction Interval Widths"
    }
  } else if (is.numeric(object)) {
    x <- as.numeric(object)
    if (is.null(coords)) {
      stop("When 'object' is a numeric vector, 'coords' must be provided.")
    }
    target_label <- "Numeric Values"
  } else {
    stop("'object' must be an object of class 'spconform' or a numeric vector.")
  }
  
  n <- length(x)
  if (n < 8) {
    stop("Anisotropy test requires at least n = 8 observations.")
  }
  
  coords <- as.matrix(coords)
  valid <- !is.na(x) & !is.na(coords[, 1]) & !is.na(coords[, 2])
  if (sum(valid) < n) {
    x <- x[valid]
    coords <- coords[valid, , drop = FALSE]
    n <- length(x)
  }
  
  # Pairwise differences and angles
  # All unique pairs (i, j) with i < j
  pair_idx <- which(lower.tri(matrix(0, n, n)), arr.ind = TRUE)
  i_idx <- pair_idx[, 1]
  j_idx <- pair_idx[, 2]
  
  dx <- coords[j_idx, 1] - coords[i_idx, 1]
  dy <- coords[j_idx, 2] - coords[i_idx, 2]
  dist_pairs <- sqrt(dx^2 + dy^2)
  
  # Filter by max_dist
  if (is.null(max_dist) || max_dist <= 0) {
    max_dist <- 0.5 * max(dist_pairs)
  }
  valid_dist <- (dist_pairs > 0) & (dist_pairs <= max_dist)
  
  dx <- dx[valid_dist]
  dy <- dy[valid_dist]
  dist_pairs <- dist_pairs[valid_dist]
  sq_diff <- (x[j_idx[valid_dist]] - x[i_idx[valid_dist]])^2
  
  # Azimuth angle in degrees: [0, 180)
  # atan2(dy, dx) returns angle in radians [-pi, pi]
  angles_rad <- atan2(dy, dx)
  angles_deg <- (angles_rad * 180 / pi) %% 180
  
  # Directional bins:
  # If n_directions == 4:
  # 0 deg (E-W): [0, 22.5) U [157.5, 180)
  # 45 deg (NE-SW): [22.5, 67.5)
  # 90 deg (N-S): [67.5, 112.5)
  # 135 deg (NW-SE): [112.5, 157.5)
  dir_labels <- c("0 deg (E-W)", "45 deg (NE-SW)", "90 deg (N-S)", "135 deg (NW-SE)")
  dir_centers <- c(0, 45, 90, 135)
  
  dir_bin <- integer(length(angles_deg))
  for (k in seq_along(angles_deg)) {
    ang <- angles_deg[k]
    if (ang < 22.5 || ang >= 157.5) {
      dir_bin[k] <- 1 # 0 deg
    } else if (ang >= 22.5 && ang < 67.5) {
      dir_bin[k] <- 2 # 45 deg
    } else if (ang >= 67.5 && ang < 112.5) {
      dir_bin[k] <- 3 # 90 deg
    } else {
      dir_bin[k] <- 4 # 135 deg
    }
  }
  dir_factor <- factor(dir_bin, levels = 1:4, labels = dir_labels)
  
  df_pairs <- data.frame(gamma = 0.5 * sq_diff, dist = dist_pairs, direction = dir_factor)
  
  # Summary table of directional semivariances
  dir_summary <- data.frame(
    Direction = dir_labels,
    Angle_Deg = dir_centers,
    Pairs = as.numeric(table(dir_factor)),
    Mean_Distance = round(as.numeric(tapply(df_pairs$dist, df_pairs$direction, mean, default = 0)), 2),
    Semivariance = round(as.numeric(tapply(df_pairs$gamma, df_pairs$direction, mean, default = 0)), 4)
  )
  
  valid_dirs <- dir_summary$Semivariance[dir_summary$Pairs > 5]
  if (length(valid_dirs) >= 2) {
    aniso_ratio <- max(valid_dirs) / (min(valid_dirs) + 1e-12)
    max_idx <- which.max(dir_summary$Semivariance)
    princ_angle <- dir_summary$Angle_Deg[max_idx]
  } else {
    aniso_ratio <- 1.0
    princ_angle <- 0
  }
  
  # ANOVA Model for Anisotropy Test:
  # Isotropic model: gamma ~ dist
  # Anisotropic model: gamma ~ dist + direction
  m_iso <- stats::lm(gamma ~ dist, data = df_pairs)
  m_aniso <- stats::lm(gamma ~ dist + direction, data = df_pairs)
  
  anv <- stats::anova(m_iso, m_aniso)
  f_stat <- anv$F[2]
  df_num <- anv$Df[2]
  df_den <- anv$Res.Df[2]
  p_asymp <- anv[["Pr(>F)"]][2]
  if (is.na(f_stat)) {
    f_stat <- 0
    p_asymp <- 1.0
  }
  
  # Permutation test: rotate/shuffle coordinate pairs
  p_perm <- NA_real_
  if (n_perm > 0) {
    perm_f <- numeric(n_perm)
    for (b in seq_len(n_perm)) {
      df_b <- df_pairs
      df_b$direction <- sample(df_pairs$direction)
      m_b <- stats::lm(gamma ~ dist + direction, data = df_b)
      anv_b <- stats::anova(m_iso, m_b)
      perm_f[b] <- max(0, anv_b$F[2], na.rm = TRUE)
    }
    p_perm <- (1 + sum(perm_f >= f_stat)) / (n_perm + 1)
  }
  
  final_p <- if (!is.na(p_perm)) p_perm else p_asymp
  
  decision <- if (final_p >= 0.05) {
    "Spatial isotropy confirmed (Fail to reject H0; no significant directional bias)."
  } else {
    sprintf("Directional anisotropy detected (Reject H0; principal orientation at %d deg).", princ_angle)
  }
  
  res <- list(
    F_statistic = round(f_stat, 4),
    df = c(df_num, df_den),
    p_value = final_p,
    p_asymptotic = p_asymp,
    p_perm = p_perm,
    anisotropy_ratio = round(aniso_ratio, 3),
    principal_angle = princ_angle,
    directional_table = dir_summary,
    decision = decision,
    target = target_label,
    n_pairs = nrow(df_pairs),
    max_dist = max_dist
  )
  class(res) <- "spconform_anisotropy"
  return(res)
}

#' @export
print.spconform_anisotropy <- function(x, ...) {
  cat("-----------------------------------------------------------\n")
  cat("  spconform: Directional Anisotropy Test\n")
  cat("-----------------------------------------------------------\n")
  cat("Diagnostic Target:            ", x$target, "\n", sep = "")
  cat(sprintf("Number of Coordinate Pairs:   %d (max dist = %.2f)\n", x$n_pairs, x$max_dist))
  cat(sprintf("Test Statistic (ANOVA F):     %.4f (df1 = %d, df2 = %d)\n", 
              x$F_statistic, x$df[1], x$df[2]))
  cat("p-value:                      ", format.pval(x$p_value, eps = 1e-4, digits = 4), 
      if (!is.na(x$p_perm)) " (permutation exact)" else " (asymptotic)", "\n", sep = "")
  if (!is.na(x$p_perm)) {
    cat("  - Asymptotic p-value:        ", format.pval(x$p_asymptotic, eps = 1e-4, digits = 4), "\n", sep = "")
    cat("  - Permutation p-value:       ", format.pval(x$p_perm, eps = 1e-4, digits = 4), "\n", sep = "")
  }
  cat(sprintf("Geometric Anisotropy Ratio:   %.3f\n", x$anisotropy_ratio))
  cat(sprintf("Principal Continuity Axis:    %d degrees\n", x$principal_angle))
  cat("Status:                       ", if (x$p_value >= 0.05) "PASS (Spatial Isotropy holds)" else "FAIL (Directional Anisotropy detected)", "\n", sep = "")
  cat("Conclusion:                   ", x$decision, "\n", sep = "")
  cat("-----------------------------------------------------------\n")
  cat("Directional Semivariance Summary:\n")
  print(x$directional_table, row.names = FALSE)
  cat("-----------------------------------------------------------\n")
  invisible(x)
}

#' @export
plot.spconform_anisotropy <- function(x, ...) {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  
  graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  
  tab <- x$directional_table
  cols <- c("#2B83BA", "#ABDDA4", "#FDAE61", "#D7191C")
  
  # Panel 1: Directional Semivariance Barplot
  bp <- graphics::barplot(tab$Semivariance, names.arg = c("0 deg", "45 deg", "90 deg", "135 deg"),
                          col = cols, border = "gray30",
                          xlab = "Azimuth Direction", ylab = "Mean Semivariance gamma(h)",
                          main = "Directional Semivariances",
                          ylim = c(0, max(tab$Semivariance) * 1.25))
  graphics::grid(nx = NA, ny = NULL, col = "gray85")
  graphics::text(bp, tab$Semivariance, labels = round(tab$Semivariance, 3), pos = 3, cex = 0.9)
  
  # Panel 2: Polar Radar / Rose Diagram of Anisotropy
  angles_rad <- c(tab$Angle_Deg * pi / 180, (tab$Angle_Deg[1] + 180) * pi / 180)
  radii <- c(tab$Semivariance, tab$Semivariance[1])
  
  # Convert to Cartesian for simple base R polar plotting
  xs <- radii * cos(angles_rad)
  ys <- radii * sin(angles_rad)
  
  lim <- max(radii) * 1.3
  plot(0, 0, type = "n", xlim = c(-lim, lim), ylim = c(-lim, lim),
       xlab = "East - West", ylab = "North - South",
       main = "Anisotropy Polar Ellipse", asp = 1, axes = FALSE)
  
  # Circular grid lines
  for (r_grid in seq(lim / 3, lim, length.out = 3)) {
    ang_seq <- seq(0, 2 * pi, length.out = 100)
    graphics::lines(r_grid * cos(ang_seq), r_grid * sin(ang_seq), lty = 3, col = "gray80")
    graphics::text(r_grid * 0.7, 0, labels = round(r_grid, 2), cex = 0.7, col = "gray50")
  }
  graphics::abline(h = 0, v = 0, lty = 2, col = "gray60")
  
  # Draw polygon
  all_xs <- c(xs, -xs)
  all_ys <- c(ys, -ys)
  graphics::polygon(all_xs, all_ys, col = grDevices::adjustcolor("#2B83BA", alpha.f = 0.35),
                    border = "#2B83BA", lwd = 2)
  graphics::points(all_xs, all_ys, pch = 19, col = "#D7191C", cex = 1.2)
  
  graphics::text(lim * 0.9, 0, "East (0 deg)", pos = 4, cex = 0.8)
  graphics::text(0, lim * 0.9, "North (90 deg)", pos = 3, cex = 0.8)
}
