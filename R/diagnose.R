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
#' @param plot Logical; if \code{TRUE} (default), generates a multi-panel diagnostic plot.
#' @param ... Additional graphical parameters passed to internal plotting methods.
#'
#' @details
#' The diagnostic suite performs five rigorous audits on the prediction intervals:
#' \enumerate{
#'   \item \strong{Marginal Validity & Sharpness}: Evaluates empirical coverage
#'     \eqn{\frac{1}{m}\sum \mathbf{1}(Y_i \in C(s_i))}, mean/median interval width,
#'     and the strictly proper \strong{Winkler Interval Score (WIS)}:
#'     \deqn{\text{WIS}_\alpha(l, u, y) = (u - l) + \frac{2}{\alpha}(l - y)\mathbf{1}(y < l) + \frac{2}{\alpha}(y - u)\mathbf{1}(y > u)}
#'   \item \strong{Conditional Spatial Stratification}: Partitions the 2D spatial domain
#'     into \eqn{\text{n\_bins} \times \text{n\_bins}} equal-area quadrants and assesses local coverage.
#'   \item \strong{Convex Hull Boundary Proximity}: Assesses whether edge-effect
#'     extrapolations degrade coverage near the boundary of the spatial domain.
#'   \item \strong{Spatial Autocorrelation (Moran's I)}: Evaluates whether prediction
#'     miscoverages or nonconformity scores exhibit residual spatial clustering via Moran's \eqn{I}.
#'   \item \strong{Score Distribution}: Analyzes empirical quantiles and normality/QQ structure.
#' }
#'
#' @return An S3 object of class \code{"spconform_diagnose"} containing:
#' \describe{
#'   \item{marginal}{List with coverage, mean width, median width, sd width, Winkler score, and sample size.}
#'   \item{conditional}{Data frame of coverage and widths partitioned by spatial quadrant.}
#'   \item{boundary}{List comparing empirical metrics between interior and boundary units.}
#'   \item{spatial_autocorr}{List containing Moran's I statistic, expected value, z-score, and p-value.}
#'   \item{scores}{Summary of nonconformity score moments and quantiles.}
#'   \item{alpha}{Miscoverage level of the evaluated object.}
#' }
#'
#' @references
#' Winkler, R. L. (1972). "A Decision-Theoretic Approach to Interval Estimation."
#' \emph{Journal of the American Statistical Association}, 67(337), 187-191.
#'
#' Mao, R., Martin, R., and Reich, B. J. (2023). "Valid Conformal Prediction
#' for Dependent Data." \emph{Journal of the American Statistical Association},
#' \doi{10.1080/01621459.2022.2147531}.
#'
#' @seealso \code{\link{scp_geostatistical}}, \code{\link{scp_areal}}, \code{\link{coverage_report}}
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
#' @importFrom stats median sd quantile qqnorm qqline pnorm dist
#' @importFrom graphics abline barplot layout par text legend
#' @importFrom grDevices chull
#' @export
diagnose <- function(object, y_true, s_test = NULL, n_bins = 4, plot = TRUE, ...) {

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

  # Winkler Interval Score (WIS)
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

  ## ---- 2. Conditional Coverage by Spatial Strata ----
  conditional <- NULL
  if (!is.null(s_test) && ncol(s_test) >= 2) {
    x_q <- cut(s_test[, 1], breaks = n_bins, include.lowest = TRUE, labels = FALSE)
    y_q <- cut(s_test[, 2], breaks = n_bins, include.lowest = TRUE, labels = FALSE)
    bin_id <- paste0("Q", x_q, "-", y_q)

    conditional <- data.frame(
      bin          = character(),
      n            = integer(),
      coverage     = numeric(),
      mean_width   = numeric(),
      mean_winkler = numeric(),
      stringsAsFactors = FALSE
    )

    for (b in sort(unique(bin_id))) {
      idx <- which(bin_id == b)
      conditional <- rbind(conditional, data.frame(
        bin          = b,
        n            = length(idx),
        coverage     = mean(covered[idx], na.rm = TRUE),
        mean_width   = mean(widths[idx], na.rm = TRUE),
        mean_winkler = mean(wis_vec[idx], na.rm = TRUE)
      ))
    }
  }

  ## ---- 3. Boundary Effects (Convex Hull) ----
  boundary <- NULL
  if (!is.null(s_test) && ncol(s_test) >= 2 && identical(object$type, "geostatistical")) {
    hull <- grDevices::chull(s_test[, 1], s_test[, 2])
    hull_pts <- s_test[hull, , drop = FALSE]

    dist_to_boundary <- apply(s_test, 1, function(pt) {
      min(sqrt((hull_pts[, 1] - pt[1])^2 + (hull_pts[, 2] - pt[2])^2))
    })

    boundary_thresh <- stats::median(dist_to_boundary, na.rm = TRUE)
    is_near <- dist_to_boundary <= boundary_thresh

    boundary <- list(
      threshold     = boundary_thresh,
      near_boundary = list(
        n            = sum(is_near),
        coverage     = mean(covered[is_near], na.rm = TRUE),
        mean_width   = mean(widths[is_near], na.rm = TRUE),
        mean_winkler = mean(wis_vec[is_near], na.rm = TRUE)
      ),
      far_boundary = list(
        n            = sum(!is_near),
        coverage     = mean(covered[!is_near], na.rm = TRUE),
        mean_width   = mean(widths[!is_near], na.rm = TRUE),
        mean_winkler = mean(wis_vec[!is_near], na.rm = TRUE)
      )
    )
  }

  ## ---- 4. Spatial Autocorrelation (Moran's I on Residuals) ----
  spatial_autocorr <- NULL
  if (!is.null(s_test) && n >= 10 && ncol(s_test) >= 2) {
    raw_resids <- y_true - object$pred
    z_res <- raw_resids - mean(raw_resids)
    
    # Inverse distance spatial weight matrix
    d_mat <- as.matrix(stats::dist(s_test))
    diag(d_mat) <- Inf
    W <- 1 / d_mat
    W[is.infinite(W)] <- 0
    # Row normalize
    row_sums <- rowSums(W)
    row_sums[row_sums == 0] <- 1
    W_norm <- W / row_sums
    
    s0_w <- sum(W_norm)
    numerator <- sum(W_norm * outer(z_res, z_res))
    denominator <- sum(z_res^2)
    
    if (denominator > 0 && s0_w > 0) {
      moran_i <- (n / s0_w) * (numerator / denominator)
      exp_i <- -1 / (n - 1)
      # Approximate variance under randomization
      var_i <- (n^2) / (s0_w^2 * (n^2 - 1))
      z_score <- (moran_i - exp_i) / sqrt(pmax(1e-8, var_i))
      p_val <- 2 * (1 - stats::pnorm(abs(z_score)))
      
      spatial_autocorr <- list(
        moran_i   = moran_i,
        expected  = exp_i,
        z_score   = z_score,
        p_value   = p_val,
        is_random = (p_val >= 0.05)
      )
    }
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

  ## ---- Return S3 Object ----
  out <- list(
    marginal         = marginal,
    conditional      = conditional,
    boundary         = boundary,
    spatial_autocorr = spatial_autocorr,
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
  conditional <- x$conditional
  boundary <- x$boundary
  scores <- x$raw_scores

  n_panels <- 2
  if (!is.null(conditional)) n_panels <- n_panels + 1
  if (!is.null(boundary)) n_panels <- n_panels + 1

  layout_mat <- matrix(1:n_panels, ncol = 2, byrow = TRUE)
  if (n_panels %% 2 == 1) {
    layout_mat <- rbind(layout_mat, c(n_panels, n_panels))
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
  if (!is.null(conditional)) {
    graphics::barplot(conditional$coverage, names.arg = conditional$bin,
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


#' @export
print.spconform_diagnose <- function(x, ...) {
  nominal <- 1 - x$alpha
  cov_val <- x$marginal$coverage

  cat("======================================================================\n")
  cat("             spconform Comprehensive Diagnostic Audit Report          \n")
  cat("======================================================================\n\n")

  # 1. Marginal Status
  status_badge <- if (cov_val >= nominal) "[PASS]" else "[WARN]"
  cat(">> 1. Marginal Validity & Prediction Sharpness:\n")
  cat(sprintf("   * %s Empirical Coverage : %.3f (Target Nominal >= %.3f)\n",
              status_badge, cov_val, nominal))
  cat(sprintf("   * Mean Interval Width    : %.4f (Median = %.4f, SD = %.4f)\n",
              x$marginal$mean_width, x$marginal$median_width, x$marginal$sd_width))
  cat(sprintf("   * Winkler Interval Score : %.4f (Strictly Proper Loss)\n",
              x$marginal$mean_winkler))
  cat(sprintf("   * Total Target Units     : n = %d (Covered = %d, Miscovered = %d)\n\n",
              x$marginal$n, x$marginal$n_covered, x$marginal$n - x$marginal$n_covered))

  # 2. Spatial Autocorrelation
  if (!is.null(x$spatial_autocorr)) {
    auto <- x$spatial_autocorr
    auto_status <- if (auto$is_random) "[PASS]" else "[NOTE]"
    cat(">> 2. Spatial Residual Autocorrelation (Moran's I Audit):\n")
    cat(sprintf("   * %s Moran's I Statistic: %.4f (Expected = %.4f, z = %.2f, p-value = %.4f)\n",
                auto_status, auto$moran_i, auto$expected, auto$z_score, auto$p_value))
    if (auto$is_random) {
      cat("   * Conclusion: Residuals show NO significant spatial clustering (Independent).\n\n")
    } else {
      cat("   * Conclusion: Moderate spatial residual structure detected; localized weights active.\n\n")
    }
  }

  # 3. Spatial Strata Quadrants
  if (!is.null(x$conditional)) {
    cat(">> 3. Conditional Coverage across Spatial Quadrants:\n")
    for (i in seq_len(nrow(x$conditional))) {
      cat(sprintf("   - Strata %-7s: Cov = %5.1f%% | Mean Width = %6.3f | WIS = %6.3f (n = %d)\n",
                  x$conditional$bin[i],
                  x$conditional$coverage[i] * 100,
                  x$conditional$mean_width[i],
                  x$conditional$mean_winkler[i],
                  x$conditional$n[i]))
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
  cat(sprintf("   * Mean = %.4f | Median = %.4f | SD = %.4f | Q90 = %.4f | Q95 = %.4f\n",
              x$scores$mean, x$scores$median, x$scores$sd, x$scores$q90, x$scores$q95))
  cat("======================================================================\n")

  invisible(x)
}