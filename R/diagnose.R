#' Spatial Autocorrelation Diagnostics (Moran's I) for Conformal Residuals
#'
#' Evaluates spatial autocorrelation and clustering among prediction residuals
#' and binary miscoverage indicators using Moran's \eqn{I} statistic.
#'
#' @param object An object of class \code{"spconform"} or numeric vector of predictions.
#' @param y_true Numeric vector of true observed responses.
#' @param s_test Numeric matrix of spatial coordinates (\eqn{m \times p}).
#'
#' @return An S3 object of class \code{"spconform_moran"} containing:
#' \describe{
#'   \item{moran_i}{Observed Moran's I statistic for prediction residuals.}
#'   \item{expected}{Theoretical expectation \eqn{-1/(n-1)} under spatial independence.}
#'   \item{z_score}{Standardized normal z-score under the randomization assumption.}
#'   \item{p_value}{Two-sided p-value against the null hypothesis of spatial randomness.}
#'   \item{is_random}{Logical; \code{TRUE} if p-value >= 0.05.}
#'   \item{moran_i_miscoverage}{Moran's I on binary hit/miss indicators if available.}
#' }
#'
#' @references
#' Cliff, A. D. and Ord, J. K. (1981). \emph{Spatial Processes: Models & Applications}. Pion, London.
#'
#' Mao, H., Martin, R., and Reich, B. J. (2024). Valid Conformal Prediction
#' for Dependent Data. \emph{Journal of the American Statistical Association},
#' 119(546), 1269--1281. \doi{10.1080/01621459.2022.2147531}.
#'
#' @examples
#' set.seed(42)
#' n <- 50
#' s <- matrix(runif(n * 2), n, 2)
#' y <- sin(3 * s[, 1]) + rnorm(n, sd = 0.2)
#' pfun <- function(s_tr, y_tr, s_new) rep(mean(y_tr), nrow(s_new))
#' fit <- scp_geostatistical(s[1:35, ], y[1:35], s[36:50, ], pfun, alpha = 0.1)
#' m_res <- spatial_moran(fit, y[36:50], s[36:50, ])
#' print(m_res)
#'
#' @importFrom stats dist pnorm
#' @export
spatial_moran <- function(object, y_true, s_test) {
  if (inherits(object, "spconform")) {
    pred <- as.vector(object$pred)
    covered <- as.vector((as.vector(y_true) >= as.vector(object$lower)) & (as.vector(y_true) <= as.vector(object$upper)))
  } else {
    pred <- as.vector(as.numeric(object))
    covered <- NULL
  }

  y_true <- as.vector(as.numeric(y_true))
  s_test <- as.matrix(s_test)
  n <- length(y_true)

  if (length(pred) != n) {
    stop("length of predictions must equal length(y_true)")
  }
  if (nrow(s_test) != n) {
    stop("nrow(s_test) must equal length(y_true)")
  }
  if (n < 5) {
    stop("Need at least 5 spatial locations to compute Moran's I.")
  }

  raw_resids <- as.vector(y_true - pred)
  z_res <- as.vector(raw_resids - mean(raw_resids))

  # Inverse distance spatial weight matrix
  d_mat <- as.matrix(stats::dist(s_test))
  diag(d_mat) <- Inf
  W <- 1 / d_mat
  W[is.infinite(W)] <- 0
  row_sums <- rowSums(W)
  row_sums[row_sums == 0] <- 1
  W_norm <- W / row_sums

  s0_w <- sum(W_norm)
  numerator <- sum(W_norm * outer(z_res, z_res))
  denominator <- sum(z_res^2)

  moran_i <- 0
  exp_i <- -1 / (n - 1)
  z_score <- 0
  p_val <- 1

  if (denominator > 0 && s0_w > 0) {
    moran_i <- (n / s0_w) * (numerator / denominator)
    var_i <- (n^2) / (s0_w^2 * (n^2 - 1))
    z_score <- (moran_i - exp_i) / sqrt(pmax(1e-8, var_i))
    p_val <- 2 * (1 - stats::pnorm(abs(z_score)))
  }

  # Moran on binary miscoverage hits
  moran_hit <- NULL
  if (!is.null(covered)) {
    z_hit <- as.vector(as.numeric(covered) - mean(covered))
    denom_hit <- sum(z_hit^2)
    if (denom_hit > 0 && s0_w > 0) {
      moran_hit <- (n / s0_w) * (sum(W_norm * outer(z_hit, z_hit)) / denom_hit)
    }
  }

  out <- list(
    moran_i              = moran_i,
    expected             = exp_i,
    z_score              = z_score,
    p_value              = p_val,
    is_random            = (p_val >= 0.05),
    moran_i_miscoverage  = moran_hit,
    n                    = n
  )
  class(out) <- "spconform_moran"
  out
}

#' Print Method for spconform_moran Objects
#'
#' @param x An object of class \code{"spconform_moran"}.
#' @param ... Further arguments (unused).
#' @export
print.spconform_moran <- function(x, ...) {
  cat("Spatial Autocorrelation Audit (Moran's I on Residuals)\n")
  cat("----------------------------------------------------\n")
  cat(sprintf("Observed Moran's I : %+.4f\n", x$moran_i))
  cat(sprintf("Expected Value      : %+.4f\n", x$expected))
  cat(sprintf("Standard Normal Z   : %+.4f\n", x$z_score))
  cat(sprintf("Two-sided p-value   : %.4f %s\n", x$p_value,
              if (x$is_random) "(No residual spatial clustering, p >= 0.05)"
              else "(Significant spatial clustering detected, p < 0.05)"))
  if (!is.null(x$moran_i_miscoverage)) {
    cat(sprintf("Moran's I on Hits   : %+.4f\n", x$moran_i_miscoverage))
  }
  invisible(x)
}

#' Domain Boundary Proximity Diagnostics for Conformal Prediction
#'
#' Evaluates whether prediction intervals degrade near the convex-hull boundary
#' of the spatial prediction domain versus interior locations.
#'
#' @param object An object of class \code{"spconform"}.
#' @param y_true Numeric vector of true observed responses.
#' @param s_test Numeric matrix of spatial coordinates (\eqn{m \times 2}).
#'
#' @return An S3 object of class \code{"spconform_boundary"} comparing empirical coverage,
#'   interval width, and Winkler scores across edge and interior subdomains.
#'
#' @references
#' Winkler, R. L. (1972). A Decision-Theoretic Approach to Interval Estimation.
#' \emph{Journal of the American Statistical Association}, 67(337), 187--191.
#'
#' @examples
#' set.seed(42)
#' n <- 50
#' s <- matrix(runif(n * 2), n, 2)
#' y <- sin(3 * s[, 1]) + rnorm(n, sd = 0.2)
#' pfun <- function(s_tr, y_tr, s_new) rep(mean(y_tr), nrow(s_new))
#' fit <- scp_geostatistical(s[1:35, ], y[1:35], s[36:50, ], pfun, alpha = 0.1)
#' b_res <- boundary_diagnostics(fit, y[36:50], s[36:50, ])
#' print(b_res)
#'
#' @importFrom grDevices chull
#' @importFrom stats median
#' @export
boundary_diagnostics <- function(object, y_true, s_test) {
  if (!inherits(object, "spconform")) stop("'object' must be of class 'spconform'")
  s_test <- as.matrix(s_test)
  y_true <- as.numeric(y_true)
  n <- length(object$pred)

  if (ncol(s_test) < 2) stop("'s_test' must have at least 2 coordinate columns.")
  if (nrow(s_test) != n) stop("nrow(s_test) must equal length(object$pred)")

  alpha <- object$alpha
  widths <- object$upper - object$lower
  covered <- (y_true >= object$lower) & (y_true <= object$upper)
  underage <- pmax(0, object$lower - y_true)
  overage  <- pmax(0, y_true - object$upper)
  wis_vec  <- widths + (2 / alpha) * underage + (2 / alpha) * overage

  hull <- grDevices::chull(s_test[, 1], s_test[, 2])
  hull_pts <- s_test[hull, , drop = FALSE]

  dist_to_boundary <- apply(s_test, 1, function(pt) {
    min(sqrt((hull_pts[, 1] - pt[1])^2 + (hull_pts[, 2] - pt[2])^2))
  })

  boundary_thresh <- stats::median(dist_to_boundary, na.rm = TRUE)
  is_near <- dist_to_boundary <= boundary_thresh

  out <- list(
    threshold     = boundary_thresh,
    near_boundary = list(
      n            = sum(is_near),
      coverage     = mean(covered[is_near], na.rm = TRUE),
      mean_width   = mean(widths[is_near], na.rm = TRUE),
      mean_winkler = mean(wis_vec[is_near], na.rm = TRUE)
    ),
    far_boundary  = list(
      n            = sum(!is_near),
      coverage     = mean(covered[!is_near], na.rm = TRUE),
      mean_width   = mean(widths[!is_near], na.rm = TRUE),
      mean_winkler = mean(wis_vec[!is_near], na.rm = TRUE)
    ),
    alpha         = alpha
  )
  class(out) <- "spconform_boundary"
  out
}

#' Print Method for spconform_boundary Objects
#'
#' @param x An object of class \code{"spconform_boundary"}.
#' @param ... Further arguments (unused).
#' @export
print.spconform_boundary <- function(x, ...) {
  cat("Domain Boundary Proximity Audit (Convex Hull)\n")
  cat("---------------------------------------------\n")
  cat(sprintf("Distance Median Threshold: %.4f\n", x$threshold))
  cat(sprintf("Near Boundary (Edge) : n = %d | Cov = %5.1f%% | Mean Width = %6.3f | WIS = %6.3f\n",
              x$near_boundary$n, x$near_boundary$coverage * 100,
              x$near_boundary$mean_width, x$near_boundary$mean_winkler))
  cat(sprintf("Far Boundary (Core) : n = %d | Cov = %5.1f%% | Mean Width = %6.3f | WIS = %6.3f\n",
              x$far_boundary$n, x$far_boundary$coverage * 100,
              x$far_boundary$mean_width, x$far_boundary$mean_winkler))
  invisible(x)
}

#' Plot Method for spconform_boundary Objects
#'
#' @param x An object of class \code{"spconform_boundary"}.
#' @param ... Further graphical arguments.
#' @importFrom graphics barplot abline
#' @export
plot.spconform_boundary <- function(x, ...) {
  nominal <- 1 - x$alpha
  graphics::barplot(c(x$near_boundary$coverage, x$far_boundary$coverage),
                    names.arg = c("Near Boundary", "Far (Interior)"),
                    col = c("coral", "steelblue"), ylab = "Coverage Rate",
                    main = "Coverage vs Boundary Proximity",
                    ylim = c(0, 1.15), ...)
  graphics::abline(h = nominal, col = "red", lty = 2, lwd = 1.5)
  invisible(x)
}

#' Conditional Spatial Strata Diagnostics for Conformal Prediction
#'
#' Partitions the spatial domain into coordinate grid strata/quadrants and audits
#' local empirical coverage, interval widths, and Winkler scores.
#'
#' @param object An object of class \code{"spconform"}.
#' @param y_true Numeric vector of true observed responses.
#' @param s_test Numeric matrix of spatial coordinates (\eqn{m \times 2}).
#' @param n_bins Integer; number of spatial bins per coordinate axis (default 4).
#'
#' @return An S3 object of class \code{"spconform_strata"} containing a data frame
#'   of strata evaluations and summary metrics.
#'
#' @references
#' Winkler, R. L. (1972). A Decision-Theoretic Approach to Interval Estimation.
#' \emph{Journal of the American Statistical Association}, 67(337), 187--191.
#'
#' @examples
#' set.seed(42)
#' n <- 50
#' s <- matrix(runif(n * 2), n, 2)
#' y <- sin(3 * s[, 1]) + rnorm(n, sd = 0.2)
#' pfun <- function(s_tr, y_tr, s_new) rep(mean(y_tr), nrow(s_new))
#' fit <- scp_geostatistical(s[1:35, ], y[1:35], s[36:50, ], pfun, alpha = 0.1)
#' s_res <- strata_diagnostics(fit, y[36:50], s[36:50, ], n_bins = 2)
#' print(s_res)
#'
#' @export
strata_diagnostics <- function(object, y_true, s_test, n_bins = 4) {
  if (!inherits(object, "spconform")) stop("'object' must be of class 'spconform'")
  s_test <- as.matrix(s_test)
  y_true <- as.numeric(y_true)
  n <- length(object$pred)

  if (ncol(s_test) < 2) stop("'s_test' must have at least 2 coordinate columns.")
  if (nrow(s_test) != n) stop("nrow(s_test) must equal length(object$pred)")

  alpha <- object$alpha
  widths <- object$upper - object$lower
  covered <- (y_true >= object$lower) & (y_true <= object$upper)
  underage <- pmax(0, object$lower - y_true)
  overage  <- pmax(0, y_true - object$upper)
  wis_vec  <- widths + (2 / alpha) * underage + (2 / alpha) * overage

  x_q <- cut(s_test[, 1], breaks = n_bins, include.lowest = TRUE, labels = FALSE)
  y_q <- cut(s_test[, 2], breaks = n_bins, include.lowest = TRUE, labels = FALSE)
  bin_id <- paste0("Q", x_q, "-", y_q)

  df <- data.frame(
    bin          = character(),
    n            = integer(),
    coverage     = numeric(),
    mean_width   = numeric(),
    mean_winkler = numeric(),
    stringsAsFactors = FALSE
  )

  for (b in sort(unique(bin_id))) {
    idx <- which(bin_id == b)
    df <- rbind(df, data.frame(
      bin          = b,
      n            = length(idx),
      coverage     = mean(covered[idx], na.rm = TRUE),
      mean_width   = mean(widths[idx], na.rm = TRUE),
      mean_winkler = mean(wis_vec[idx], na.rm = TRUE),
      stringsAsFactors = FALSE
    ))
  }

  out <- list(
    strata = df,
    alpha  = alpha,
    n_bins = n_bins
  )
  class(out) <- "spconform_strata"
  out
}

#' Print Method for spconform_strata Objects
#'
#' @param x An object of class \code{"spconform_strata"}.
#' @param ... Further arguments (unused).
#' @export
print.spconform_strata <- function(x, ...) {
  cat("Conditional Spatial Strata Audit\n")
  cat("--------------------------------\n")
  for (i in seq_len(nrow(x$strata))) {
    cat(sprintf("   - Strata %-7s: Cov = %5.1f%% | Mean Width = %6.3f | WIS = %6.3f (n = %d)\n",
                x$strata$bin[i],
                x$strata$coverage[i] * 100,
                x$strata$mean_width[i],
                x$strata$mean_winkler[i],
                x$strata$n[i]))
  }
  invisible(x)
}

#' Plot Method for spconform_strata Objects
#'
#' @param x An object of class \code{"spconform_strata"}.
#' @param ... Further graphical arguments.
#' @importFrom graphics barplot abline
#' @export
plot.spconform_strata <- function(x, ...) {
  nominal <- 1 - x$alpha
  graphics::barplot(x$strata$coverage, names.arg = x$strata$bin,
                    col = "steelblue", ylab = "Coverage Rate",
                    main = "Conditional Coverage by Spatial Strata",
                    ylim = c(0, 1.15), las = 2, cex.names = 0.7, ...)
  graphics::abline(h = nominal, col = "red", lty = 2, lwd = 1.5)
  invisible(x)
}

#' Comprehensive Diagnostic Suite for Spatial Conformal Prediction Objects
#'
#' Produces a comprehensive multi-panel diagnostic report evaluating marginal
#' coverage validity, Winkler Interval Score (WIS) sharpness, conditional coverage
#' across spatial strata, boundary proximity effects, spatial residual autocorrelation
#' (Moran's I), and nonconformity score distributions.
#'
#' @param object An object of class \code{"spconform"}, typically output from
#'   \code{\link{scp_geostatistical}} or \code{\link{scp_areal}}.
#' @param y_true Numeric vector of true observed responses at target prediction locations.
#'   Must have the same length as \code{object$pred}.
#' @param s_test Optional numeric matrix of target prediction coordinates
#'   (\eqn{m \times p}). Required for spatial stratification, boundary effects,
#'   and spatial autocorrelation diagnostics.
#' @param n_bins Integer; number of spatial bins per dimension for conditional coverage (default 4).
#' @param run_advanced Logical; if \code{TRUE} (default), executes the advanced spatial tests suite.
#' @param plot Logical; if \code{TRUE} (default), generates a multi-panel diagnostic plot.
#' @param ... Additional graphical parameters passed to internal plotting methods.
#'
#' @details
#' The diagnostic suite coordinates modular audits across the prediction intervals:
#' \enumerate{
#'   \item \strong{Marginal Validity & Sharpness}: Evaluates empirical coverage
#'     and the strictly proper Winkler Interval Score (WIS; Winkler 1972).
#'   \item \strong{Conditional Spatial Stratification}: Calls \code{\link{strata_diagnostics}}
#'     to evaluate local coverage across quadrants.
#'   \item \strong{Convex Hull Boundary Proximity}: Calls \code{\link{boundary_diagnostics}}
#'     to test whether edge-effect extrapolations degrade calibration.
#'   \item \strong{Spatial Autocorrelation (Moran's I)}: Calls \code{\link{spatial_moran}}
#'     to test for residual spatial dependence.
#'   \item \strong{Score Distribution}: Analyzes empirical quantiles and normality/QQ structure.
#' }
#'
#' @return An S3 object of class \code{"spconform_diagnose"} containing:
#' \describe{
#'   \item{marginal}{List with coverage, mean width, median width, sd width, Winkler score, and sample size.}
#'   \item{conditional}{Data frame of coverage and widths partitioned by spatial quadrant (class \code{"spconform_strata"}).}
#'   \item{boundary}{Object of class \code{"spconform_boundary"} comparing edge vs interior units.}
#'   \item{spatial_autocorr}{Object of class \code{"spconform_moran"} with Moran's I statistics.}
#'   \item{scores}{Summary of nonconformity score moments and quantiles.}
#'   \item{raw_scores}{Numeric vector of nonconformity scores.}
#'   \item{alpha}{Miscoverage level of the evaluated object.}
#' }
#'
#' @references
#' Winkler, R. L. (1972). A Decision-Theoretic Approach to Interval Estimation.
#' \emph{Journal of the American Statistical Association}, 67(337), 187--191.
#'
#' Mao, H., Martin, R., and Reich, B. J. (2024). Valid Conformal Prediction
#' for Dependent Data. \emph{Journal of the American Statistical Association},
#' 119(546), 1269--1281. \doi{10.1080/01621459.2022.2147531}.
#'
#' Susmann, H. P. (2024). Adaptive Conformal Prediction for Spatial Data.
#' \emph{Computo}, \doi{10.57750/edan-5f53}.
#'
#' @seealso \code{\link{spatial_moran}}, \code{\link{boundary_diagnostics}},
#'   \code{\link{strata_diagnostics}}, \code{\link{scp_geostatistical}}, \code{\link{scp_areal}}
#'
#' @examples
#' set.seed(42)
#' n <- 100
#' s_tr <- matrix(runif(n * 2), n, 2)
#' y_tr <- sin(3 * s_tr[, 1]) + rnorm(n, sd = 0.2)
#' s_te <- matrix(runif(30 * 2), 30, 2)
#' y_te <- sin(3 * s_te[, 1]) + rnorm(30, sd = 0.2)
#'
#' pfun <- function(s_tr, y_tr, s_new) rep(mean(y_tr), nrow(s_new))
#' out <- scp_geostatistical(s_tr, y_tr, s_te, pfun, alpha = 0.1)
#'
#' diag_res <- diagnose(out, y_true = y_te, s_test = s_te, plot = FALSE)
#' print(diag_res)
#'
#' # Programmatic access to modular components
#' as.data.frame(diag_res)
#'
#' @importFrom stats median sd quantile qqnorm qqline pnorm dist
#' @importFrom graphics abline barplot layout par text legend
#' @importFrom grDevices chull
#' @export
diagnose <- function(object, y_true, s_test = NULL, n_bins = 4, run_advanced = TRUE, plot = TRUE, ...) {

  if (!inherits(object, "spconform")) {
    stop("'object' must be of class 'spconform'")
  }

  n <- length(object$pred)
  y_true <- as.numeric(y_true)
  object$pred <- as.numeric(object$pred)
  object$lower <- as.numeric(object$lower)
  object$upper <- as.numeric(object$upper)
  if (length(y_true) != n) {
    stop("length(y_true) must equal length(object$pred)")
  }

  alpha <- object$alpha
  nominal <- 1 - alpha

  ## ---- 1. Marginal Coverage & Winkler Interval Score ----
  covered <- (y_true >= object$lower) & (y_true <= object$upper)
  widths  <- object$upper - object$lower

  underage <- pmax(0, object$lower - y_true)
  overage  <- pmax(0, y_true - object$upper)
  wis_vec  <- widths + (2 / alpha) * underage + (2 / alpha) * overage

  marginal <- list(
    coverage     = mean(covered, na.rm = TRUE),
    mean_width   = mean(widths, na.rm = TRUE),
    median_width = stats::median(widths, na.rm = TRUE),
    sd_width     = stats::sd(widths, na.rm = TRUE),
    mean_winkler = mean(wis_vec, na.rm = TRUE),
    n            = n,
    n_covered    = sum(covered, na.rm = TRUE)
  )

  ## ---- 2. Modular Conditional Strata Diagnostics ----
  conditional <- NULL
  if (!is.null(s_test) && ncol(s_test) >= 2) {
    conditional <- strata_diagnostics(object, y_true, s_test, n_bins = n_bins)
  }

  ## ---- 3. Modular Boundary Proximity Diagnostics ----
  boundary <- NULL
  if (!is.null(s_test) && ncol(s_test) >= 2 && identical(object$type, "geostatistical")) {
    boundary <- boundary_diagnostics(object, y_true, s_test)
  }

  ## ---- 4. Modular Spatial Autocorrelation Diagnostics ----
  spatial_autocorr <- NULL
  if (!is.null(s_test) && n >= 5 && ncol(s_test) >= 2) {
    spatial_autocorr <- spatial_moran(object, y_true, s_test)
  }

  ## ---- 4b. Advanced Spatial Diagnostic Tests (v0.2.2) ----
  conditional_test <- NULL
  geary_res        <- NULL
  lisa_res         <- NULL
  aniso_res        <- NULL

  if (isTRUE(run_advanced) && !is.null(s_test) && n >= 10 && ncol(s_test) >= 2) {
    resids_vec <- if (inherits(object, "spconform")) as.vector(y_true - object$pred) else as.vector(y_true - object)
    conditional_test <- tryCatch(spatial_conditional_coverage(object, y_true = y_true, coords = s_test), error = function(e) NULL)
    geary_res        <- tryCatch(spatial_geary(resids_vec, coords = s_test), error = function(e) NULL)
    lisa_res         <- tryCatch(spatial_lisa(object, coords = s_test, y_true = y_true, target = "residuals"), error = function(e) NULL)
    aniso_res        <- tryCatch(spatial_anisotropy(object, coords = s_test, y_true = y_true), error = function(e) NULL)
  }

  ## ---- 5. Nonconformity Score Distribution ----
  scores <- pmax(object$upper - object$pred, object$pred - object$lower)

  scores_summary <- list(
    mean   = mean(scores, na.rm = TRUE),
    median = stats::median(scores, na.rm = TRUE),
    sd     = stats::sd(scores, na.rm = TRUE),
    q90    = stats::quantile(scores, 0.90, na.rm = TRUE),
    q95    = stats::quantile(scores, 0.95, na.rm = TRUE)
  )

  out <- list(
    marginal         = marginal,
    conditional      = conditional,
    boundary         = boundary,
    spatial_autocorr = spatial_autocorr,
    conditional_test = conditional_test,
    geary            = geary_res,
    lisa             = lisa_res,
    anisotropy       = aniso_res,
    scores           = scores_summary,
    raw_scores       = scores,
    alpha            = object$alpha
  )

  class(out) <- c("spconform_diagnose", "list")

  if (isTRUE(plot)) {
    plot(out)
  }

  out
}

#' Coerce an spconform_diagnose Object to a Data Frame
#'
#' Extracts conditional spatial strata and diagnostic statistics as a tidy
#' data frame for programmatic inspection.
#'
#' @param x An object of class \code{"spconform_diagnose"}.
#' @param row.names Optional row names.
#' @param optional Logical; passed to \code{as.data.frame}.
#' @param ... Additional arguments (unused).
#'
#' @return A tidy \code{data.frame} of spatial diagnostic strata.
#' @export
as.data.frame.spconform_diagnose <- function(x, row.names = NULL, optional = FALSE, ...) {
  if (!is.null(x$conditional)) {
    df <- x$conditional$strata
  } else {
    df <- data.frame(
      metric = c("coverage", "mean_width", "mean_winkler", "n"),
      value  = c(x$marginal$coverage, x$marginal$mean_width, x$marginal$mean_winkler, x$marginal$n),
      stringsAsFactors = FALSE
    )
  }
  df
}

#' Plot Diagnostic Audit for an spconform_diagnose Object
#'
#' Generates a multi-panel visual audit of marginal coverage, nonconformity
#' score QQ-distribution, conditional spatial strata coverage, and boundary
#' proximity effects.
#'
#' @param x An object of class \code{"spconform_diagnose"}.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return Invisibly returns the input object \code{x}.
#' @export
plot.spconform_diagnose <- function(x, ...) {
  nominal <- 1 - x$alpha
  marginal <- x$marginal
  conditional_df <- if (!is.null(x$conditional)) x$conditional$strata else NULL
  boundary <- x$boundary
  scores <- x$raw_scores

  n_panels <- 2
  if (!is.null(conditional_df)) n_panels <- n_panels + 1
  if (!is.null(boundary)) n_panels <- n_panels + 1

  layout_mat <- if (n_panels == 2) {
    matrix(1:2, nrow = 1, ncol = 2)
  } else if (n_panels == 3) {
    matrix(c(1, 2, 3, 3), nrow = 2, ncol = 2, byrow = TRUE)
  } else {
    matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE)
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    graphics::layout(1)
  })

  graphics::layout(layout_mat)
  graphics::par(mar = c(4, 4, 3, 1))

  # Panel 1: Marginal Coverage
  bp <- graphics::barplot(c(marginal$coverage, nominal),
                          names.arg = c("Empirical", "Nominal"),
                          col = c("steelblue", "gray80"),
                          ylab = "Coverage Rate",
                          main = sprintf("Marginal Coverage (n = %d, WIS = %.2f)",
                                         marginal$n, marginal$mean_winkler),
                          ylim = c(0, 1.15))
  graphics::abline(h = nominal, col = "red", lty = 2, lwd = 1.5)
  graphics::text(bp[1], min(marginal$coverage + 0.06, 1.06),
                 labels = sprintf("%.3f", marginal$coverage), cex = 1.1, font = 2)
  graphics::text(bp[2], nominal + 0.06,
                 labels = sprintf("%.3f", nominal), cex = 1.1, font = 2)

  # Panel 2: QQ Plot
  if (!is.null(scores) && length(scores) > 0) {
    stats::qqnorm(scores, main = "QQ-Plot of Nonconformity Scores",
                  pch = 20, col = "steelblue")
    stats::qqline(scores, col = "red", lty = 2, lwd = 1.5)
  }

  # Panel 3: Spatial Strata Quadrants
  if (!is.null(conditional_df)) {
    graphics::barplot(conditional_df$coverage, names.arg = conditional_df$bin,
                      col = "steelblue", ylab = "Coverage Rate",
                      main = "Conditional Coverage by Spatial Strata",
                      ylim = c(0, 1.1), las = 2, cex.names = 0.7)
    graphics::abline(h = nominal, col = "red", lty = 2, lwd = 1.5)
  }

  # Panel 4: Boundary Proximity
  if (!is.null(boundary)) {
    graphics::barplot(c(boundary$near_boundary$coverage, boundary$far_boundary$coverage),
                      names.arg = c("Near Boundary", "Far (Interior)"),
                      col = c("coral", "steelblue"), ylab = "Coverage Rate",
                      main = "Coverage vs Boundary Proximity",
                      ylim = c(0, 1.1))
    graphics::abline(h = nominal, col = "red", lty = 2, lwd = 1.5)
  }

  invisible(x)
}

#' Print Diagnostic Audit for an spconform_diagnose Object
#'
#' @param x An object of class \code{"spconform_diagnose"}.
#' @param ... Further arguments (unused).
#' @export
print.spconform_diagnose <- function(x, ...) {
  nominal <- 1 - x$alpha
  cat("======================================================================\n")
  cat("  spconform Diagnostic Audit Report\n")
  cat("======================================================================\n\n")

  # 1. Marginal
  cat(">> 1. Marginal Coverage Validity & Sharpness:\n")
  cat(sprintf("   - Target Nominal Level  : %5.1f%%\n", nominal * 100))
  cat(sprintf("   - Empirical Coverage    : %5.1f%% (%d of %d covered)\n",
              x$marginal$coverage * 100, x$marginal$n_covered, x$marginal$n))
  cat(sprintf("   - Interval Width        : Mean = %.3f | Median = %.3f | SD = %.3f\n",
              x$marginal$mean_width, x$marginal$median_width, x$marginal$sd_width))
  cat(sprintf("   - Winkler Score (WIS)   : Mean = %.3f (Proper Interval Score)\n\n",
              x$marginal$mean_winkler))

  # 2. Moran
  if (!is.null(x$spatial_autocorr)) {
    cat(">> 2. Spatial Autocorrelation Audit (Moran's I on Residuals):\n")
    cat(sprintf("   - Moran's I Statistic   : %+.4f (Expected = %+.4f)\n",
                x$spatial_autocorr$moran_i, x$spatial_autocorr$expected))
    cat(sprintf("   - Standardized Z-Score  : %+.4f (p-value = %.4f)\n",
                x$spatial_autocorr$z_score, x$spatial_autocorr$p_value))
    if (x$spatial_autocorr$is_random) {
      cat("   - Residual Pattern      : Pass (No significant spatial clustering, p >= 0.05)\n\n")
    } else {
      cat("   - Residual Pattern      : Warning (Significant spatial clustering detected, p < 0.05)\n\n")
    }
  }

  # 3. Spatial Strata Quadrants
  if (!is.null(x$conditional)) {
    cat(">> 3. Conditional Coverage across Spatial Quadrants:\n")
    cond_df <- x$conditional$strata
    for (i in seq_len(nrow(cond_df))) {
      cat(sprintf("   - Strata %-7s: Cov = %5.1f%% | Mean Width = %6.3f | WIS = %6.3f (n = %d)\n",
                  cond_df$bin[i],
                  cond_df$coverage[i] * 100,
                  cond_df$mean_width[i],
                  cond_df$mean_winkler[i],
                  cond_df$n[i]))
    }
    cat("\n")
  }

  # 4. Boundary Effects
  if (!is.null(x$boundary)) {
    b <- x$boundary
    cat(">> 4. Domain Boundary Effect (Convex Hull):\n")
    cat(sprintf("   - Near Boundary (Edge) : Cov = %5.1f%% | Mean Width = %6.3f | WIS = %6.3f (n = %d)\n",
                b$near_boundary$coverage * 100,
                b$near_boundary$mean_width,
                b$near_boundary$mean_winkler,
                b$near_boundary$n))
    cat(sprintf("   - Far Boundary (Core) : Cov = %5.1f%% | Mean Width = %6.3f | WIS = %6.3f (n = %d)\n\n",
                b$far_boundary$coverage * 100,
                b$far_boundary$mean_width,
                b$far_boundary$mean_winkler,
                b$far_boundary$n))
  }

  # 5. Nonconformity Scores
  cat(">> 5. Nonconformity Score Distribution Moments:\n")
  cat(sprintf("   * Mean = %.4f | Median = %.4f | SD = %.4f | Q90 = %.4f | Q95 = %.4f\n\n",
              x$scores$mean, x$scores$median, x$scores$sd, x$scores$q90, x$scores$q95))

  # 6. Advanced Spatial Diagnostic Tests (v0.2.2)
  if (!is.null(x$conditional_test)) {
    ct <- x$conditional_test
    cat(">> 6. Spatial Conditional Coverage Omnibus Test:\n")
    cat(sprintf("   - Omnibus Chi-Square    : %.4f (df = %d, p = %.4f)\n",
                ct$statistic, ct$df, ct$p_value))
    cat(sprintf("   - Spatial Calibration   : %s\n\n",
                if (ct$p_value >= 0.05) "Pass (No significant spatial miscoverage bias)"
                else "Warning (Significant spatial miscalibration detected)"))
  }

  if (!is.null(x$lisa)) {
    ls <- x$lisa
    counts <- table(ls$quadrant)
    hh <- if ("High-High" %in% names(counts)) as.integer(counts[["High-High"]]) else 0L
    ll <- if ("Low-Low" %in% names(counts)) as.integer(counts[["Low-Low"]]) else 0L
    hl <- if ("High-Low" %in% names(counts)) as.integer(counts[["High-Low"]]) else 0L
    lh <- if ("Low-High" %in% names(counts)) as.integer(counts[["Low-High"]]) else 0L
    is_sig <- (as.character(ls$quadrant) != "Not Significant")
    n_sig <- sum(is_sig)
    pct_sig <- 100 * mean(is_sig)
    cat(">> 7. Local Indicators of Spatial Association (LISA Residual Clusters):\n")
    cat(sprintf("   - Hotspots (High-High)  : %d | Coldspots (Low-Low): %d | Outliers (HL/LH): %d/%d\n",
                hh, ll, hl, lh))
    cat(sprintf("   - Significant Localities: %d of %d (%.1f%% at alpha = %.2f, %s)\n\n",
                n_sig, ls$n, pct_sig, ls$significance_level, ls$p_adjust_method))
  }

  if (!is.null(x$anisotropy)) {
    an <- x$anisotropy
    cat(">> 8. Directional Semivariance & Spatial Anisotropy:\n")
    cat(sprintf("   - Geometric Anisotropy  : Ratio = %.2fx (%s continuity along %d deg)\n",
                an$anisotropy_ratio,
                if (isTRUE(an$is_anisotropic)) "Directional" else "Isotropic",
                an$principal_angle))
    cat(sprintf("   - Directional ANOVA     : F = %.3f (p = %.4f)\n",
                an$F_statistic, an$p_value))
  }
  cat("======================================================================\n")

  invisible(x)
}
