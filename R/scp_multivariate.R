#' Spatially Localized Multivariate Conformal Prediction
#'
#' Constructs finite-sample valid, spatially adaptive prediction ellipsoids and
#' simultaneous bounding intervals for multivariate responses \eqn{\mathbf{Y}(s) \in \mathbb{R}^d}
#' (\eqn{d \ge 2}) under spatial heteroskedasticity and localized inter-variable correlation.
#'
#' @param s_train Numeric matrix of training spatial coordinates (\eqn{n \times p}).
#' @param y_train Numeric matrix or data frame of multivariate responses (\eqn{n \times d}),
#'   where \eqn{d \ge 2}.
#' @param s_test Numeric matrix of target test coordinates (\eqn{m \times p}).
#' @param pred_fun Function with signature \code{pred_fun(s_fit, y_fit, s_eval)} returning
#'   an \eqn{m \times d} numeric matrix of point predictions.
#' @param alpha Significance level (miscoverage budget), e.g., 0.10 for 90\% joint confidence. Default 0.10.
#' @param split Proportion of training points allocated to model fitting (default 0.50).
#' @param h Spatial kernel bandwidth. If \code{NULL} (default), estimated via median pairwise distance.
#' @param kernel Spatial kernel weighting function: \code{"gaussian"} (default), \code{"exponential"}, or \code{"bisquare"}.
#' @param lambda Non-negative regularization parameter for covariance matrix inversion (default \code{1e-4}).
#'
#' @details
#' At each target prediction location \eqn{s_0}, standard multivariate conformal prediction
#' assumes exchangeable residuals and a globally homogeneous covariance matrix. In spatial domains,
#' however, both the dispersion and cross-variable correlations vary across space.
#'
#' \code{scp_multivariate} implements localized spatial Mahalanobis conformity:
#' \enumerate{
#'   \item Computes localized spatial kernel weights \eqn{w_i(s_0)} for each calibration location.
#'   \item Estimates the spatially localized covariance matrix:
#'     \deqn{\hat{\mathbf{\Sigma}}(s_0) = \sum_{i \in \text{cal}} w_i(s_0) \mathbf{e}_i \mathbf{e}_i^\top + \lambda \mathrm{tr}(\mathbf{\Sigma}_{\text{global}}) \mathbf{I}_d}
#'   \item Evaluates localized Mahalanobis distances:
#'     \deqn{R_i(s_0) = \sqrt{\mathbf{e}_i^\top \hat{\mathbf{\Sigma}}(s_0)^{-1} \mathbf{e}_i}}
#'   \item Extracts the localized \eqn{(1 - \alpha)} weighted conformal quantile \eqn{\hat{q}(s_0)}.
#'   \item Forms the \eqn{(1 - \alpha)} spatial prediction ellipsoid:
#'     \deqn{\mathcal{E}_{1-\alpha}(s_0) = [ \mathbf{y} \in \mathbb{R}^d : (\mathbf{y} - \hat{\mathbf{Y}}(s_0))^\top \hat{\mathbf{\Sigma}}(s_0)^{-1} (\mathbf{y} - \hat{\mathbf{Y}}(s_0)) \le \hat{q}(s_0)^2 ]}
#'   \item Derives simultaneous projected coordinate-wise bounding intervals:
#'     \deqn{\mathrm{PI}_k(s_0) = \left[ \hat{Y}_k(s_0) - \hat{q}(s_0) \sqrt{\hat{\Sigma}_{kk}(s_0)}, \; \hat{Y}_k(s_0) + \hat{q}(s_0) \sqrt{\hat{\Sigma}_{kk}(s_0)} \right]}
#' }
#'
#' @return An S3 object of class \code{c("spconform_multivariate", "spconform")} containing:
#' \describe{
#'   \item{pred}{Matrix (\eqn{m \times d}) of multivariate point predictions.}
#'   \item{lower}{Matrix (\eqn{m \times d}) of projected lower interval bounds.}
#'   \item{upper}{Matrix (\eqn{m \times d}) of projected upper interval bounds.}
#'   \item{cov_matrices}{List of \eqn{m} localized \eqn{d \times d} covariance matrices.}
#'   \item{quantiles}{Numeric vector of length \eqn{m} of localized conformal radii \eqn{\hat{q}(s_0)}.}
#'   \item{ellipsoid_volumes}{Numeric vector of length \eqn{m} of ellipsoid volume metrics.}
#'   \item{alpha}{Miscoverage rate.}
#'   \item{h}{Bandwidth used.}
#'   \item{dimension}{Response dimension \eqn{d}.}
#'   \item{response_names}{Names of response variables.}
#'   \item{s_test}{Test coordinates.}
#' }
#'
#' @references
#' Barber, R. F., Candes, E. J., Ramdas, A., and Tibshirani, R. J. (2021).
#' Predictive inference with the jackknife+. \emph{The Annals of Statistics}, 49(1), 486--507.
#'
#' Lei, J., G\'Sell, M., Rinaldo, A., Tibshirani, R. J., and Wasserman, L. (2018).
#' Distribution-free predictive inference for regression. \emph{Journal of the American Statistical Association}, 113(523), 1094--1111.
#'
#' @examples
#' set.seed(42)
#' n <- 80
#' s <- matrix(runif(n * 2, 0, 10), n, 2)
#' y1 <- sin(s[, 1]) + rnorm(n, sd = 0.3)
#' y2 <- cos(s[, 2]) + 0.5 * y1 + rnorm(n, sd = 0.3)
#' y <- cbind(cadmium = y1, zinc = y2)
#'
#' s_tr <- s[1:60, ]
#' y_tr <- y[1:60, ]
#' s_te <- s[61:80, ]
#'
#' pfun <- function(st, yt, sn) {
#'   m1 <- lm(yt[, 1] ~ st[, 1] + st[, 2])
#'   m2 <- lm(yt[, 2] ~ st[, 1] + st[, 2])
#'   cbind(predict(m1, newdata = as.data.frame(sn)),
#'         predict(m2, newdata = as.data.frame(sn)))
#' }
#'
#' fit <- scp_multivariate(s_tr, y_tr, s_te, pfun, alpha = 0.10)
#' print(fit)
#'
#' @importFrom stats dist cov median
#' @importFrom grDevices adjustcolor
#' @export
scp_multivariate <- function(s_train, y_train, s_test, pred_fun,
                             alpha = 0.10, split = 0.50,
                             h = NULL, kernel = c("gaussian", "exponential", "bisquare"),
                             lambda = 1e-4) {

  kernel <- match.arg(kernel)
  s_train <- as.matrix(s_train)
  s_test  <- as.matrix(s_test)
  y_train <- as.matrix(y_train)

  n <- nrow(s_train)
  m <- nrow(s_test)
  d <- ncol(y_train)

  if (d < 2) {
    stop("scp_multivariate requires at least 2 response dimensions (d >= 2).")
  }
  if (nrow(y_train) != n) {
    stop("nrow(y_train) must equal nrow(s_train).")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be strictly between 0 and 1.")
  }

  resp_names <- colnames(y_train)
  if (is.null(resp_names)) {
    resp_names <- paste0("Y", seq_len(d))
    colnames(y_train) <- resp_names
  }

  # 1. Train / Calibration split
  n_fit <- max(d + 2, floor(split * n))
  n_cal <- n - n_fit
  if (n_cal < d + 2) {
    stop("Not enough calibration samples for multivariate covariance estimation. Increase n or decrease split.")
  }

  idx_fit <- seq_len(n_fit)
  idx_cal <- seq(n_fit + 1, n)

  s_fit <- s_train[idx_fit, , drop = FALSE]
  y_fit <- y_train[idx_fit, , drop = FALSE]

  s_cal <- s_train[idx_cal, , drop = FALSE]
  y_cal <- y_train[idx_cal, , drop = FALSE]

  # 2. Fit and predict on calibration set
  pred_cal <- pred_fun(s_fit, y_fit, s_cal)
  pred_cal <- as.matrix(pred_cal)
  if (nrow(pred_cal) != n_cal || ncol(pred_cal) != d) {
    stop("pred_fun(s_fit, y_fit, s_cal) must return a matrix of dimension n_cal x d.")
  }

  resids_cal <- y_cal - pred_cal

  # Global covariance for ridge shrinkage scale
  cov_global <- stats::cov(resids_cal)
  ridge_term <- lambda * sum(diag(cov_global)) * diag(d)

  # 3. Fit and predict on test set
  pred_test <- pred_fun(s_fit, y_fit, s_test)
  pred_test <- as.matrix(pred_test)
  colnames(pred_test) <- resp_names

  # 4. Bandwidth selection
  if (is.null(h)) {
    d_all <- as.vector(stats::dist(s_train))
    d_all <- d_all[d_all > 0]
    h <- if (length(d_all) > 0) stats::median(d_all) else 1.0
  }

  # Kernel weight function
  calc_weights <- function(dists, bw, kern) {
    u <- dists / bw
    if (kern == "gaussian") {
      w <- exp(-0.5 * u^2)
    } else if (kern == "exponential") {
      w <- exp(-u)
    } else if (kern == "bisquare") {
      w <- ifelse(u <= 1, (1 - u^2)^2, 0)
    }
    w
  }

  # Distance matrix between test and calibration
  diff_x <- outer(s_test[, 1], s_cal[, 1], "-")
  diff_y <- outer(s_test[, 2], s_cal[, 2], "-")
  dist_test_cal <- sqrt(diff_x^2 + diff_y^2)

  lower_mat <- matrix(NA_real_, nrow = m, ncol = d, dimnames = list(NULL, resp_names))
  upper_mat <- matrix(NA_real_, nrow = m, ncol = d, dimnames = list(NULL, resp_names))
  cov_list  <- vector("list", m)
  quantiles <- numeric(m)
  volumes   <- numeric(m)

  # Gamma constant for d-dimensional ellipsoid volume: pi^(d/2) / Gamma(d/2 + 1)
  vol_const <- (pi^(d / 2)) / gamma(d / 2 + 1)

  for (j in seq_len(m)) {
    dists_j <- dist_test_cal[j, ]
    w_j <- calc_weights(dists_j, h, kernel)
    if (sum(w_j) <= 1e-12) {
      w_j <- rep(1 / n_cal, n_cal)
    } else {
      w_j <- w_j / sum(w_j)
    }

    # Localized weighted covariance matrix
    # Sigma(s0) = sum_i w_i e_i e_i^T + ridge
    Sigma_j <- matrix(0, d, d)
    for (i in seq_len(n_cal)) {
      e_i <- resids_cal[i, ]
      Sigma_j <- Sigma_j + w_j[i] * tcrossprod(e_i)
    }
    Sigma_j <- Sigma_j + ridge_term

    # Invert covariance matrix safely
    Sigma_inv_j <- solve(Sigma_j)

    # Localized Mahalanobis conformity scores for all calibration points
    scores_j <- numeric(n_cal)
    for (i in seq_len(n_cal)) {
      e_i <- resids_cal[i, ]
      scores_j[i] <- sqrt(as.numeric(crossprod(e_i, Sigma_inv_j %*% e_i)))
    }

    # Weighted conformal quantile
    ord_j <- order(scores_j)
    scores_sorted <- scores_j[ord_j]
    weights_sorted <- w_j[ord_j]
    cum_w <- cumsum(weights_sorted)

    target_level <- 1 - alpha
    idx_q <- which(cum_w >= target_level)[1]
    if (is.na(idx_q)) idx_q <- n_cal
    q_val <- scores_sorted[idx_q]

    quantiles[j] <- q_val
    cov_list[[j]] <- Sigma_j

    # Ellipsoid volume: V = c_d * det(Sigma)^(1/2) * q^d
    det_sigma <- max(1e-12, det(Sigma_j))
    volumes[j] <- vol_const * sqrt(det_sigma) * (q_val^d)

    # Projected marginal bounding intervals:
    se_k <- sqrt(diag(Sigma_j))
    lower_mat[j, ] <- pred_test[j, ] - q_val * se_k
    upper_mat[j, ] <- pred_test[j, ] + q_val * se_k
  }

  out <- list(
    pred              = pred_test,
    lower             = lower_mat,
    upper             = upper_mat,
    cov_matrices      = cov_list,
    quantiles         = quantiles,
    ellipsoid_volumes = volumes,
    alpha             = alpha,
    h                 = h,
    dimension         = d,
    response_names    = resp_names,
    s_test            = s_test,
    type              = "multivariate"
  )

  class(out) <- c("spconform_multivariate", "spconform")
  out
}

#' Print Method for spconform_multivariate Objects
#' @param x An object of class \code{"spconform_multivariate"}.
#' @param ... Additional arguments.
#' @export
print.spconform_multivariate <- function(x, ...) {
  cat("=== Spatially Localized Multivariate Conformal Prediction ===\n")
  cat(sprintf("  Response dimension (d) : %d (%s)\n", x$dimension, paste(x$response_names, collapse = ", ")))
  cat(sprintf("  Target locations (m)   : %d\n", nrow(x$pred)))
  cat(sprintf("  Confidence Level       : %.1f%% (alpha = %.3f)\n", 100 * (1 - x$alpha), x$alpha))
  cat(sprintf("  Spatial Bandwidth (h)  : %.3f\n", x$h))
  cat(sprintf("  Mean Ellipsoid Radius  : %.3f\n", mean(x$quantiles)))
  cat(sprintf("  Mean Ellipsoid Volume  : %.3f\n", mean(x$ellipsoid_volumes)))
  cat("\nMarginal Interval Summary (Projected Hyper-rectangle):\n")
  for (k in seq_len(x$dimension)) {
    nm <- x$response_names[k]
    w_k <- x$upper[, k] - x$lower[, k]
    cat(sprintf("  Variable '%s': Mean Width = %.3f, Median = %.3f [min: %.3f, max: %.3f]\n",
                nm, mean(w_k), stats::median(w_k), min(w_k), max(w_k)))
  }
  invisible(x)
}

#' Plot Method for spconform_multivariate Objects
#'
#' Visualizes bivariate prediction ellipsoids or multivariate interval widths.
#'
#' @param x An object of class \code{"spconform_multivariate"}.
#' @param idx Integer vector of test indices to plot ellipsoids for (default first 6).
#' @param y_test Optional matrix of observed test responses for coverage validation.
#' @param ... Additional graphical arguments.
#' @importFrom graphics par plot lines points text legend polygon
#' @importFrom grDevices adjustcolor
#' @export
plot.spconform_multivariate <- function(x, idx = NULL, y_test = NULL, ...) {
  d <- x$dimension
  m <- nrow(x$pred)
  if (is.null(idx)) {
    idx <- seq_len(min(6, m))
  }

  if (d == 2) {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))
    graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

    # Panel 1: Bivariate Ellipsoids for selected points
    all_x <- c(x$lower[idx, 1], x$upper[idx, 1])
    all_y <- c(x$lower[idx, 2], x$upper[idx, 2])
    if (!is.null(y_test)) {
      all_x <- c(all_x, y_test[idx, 1])
      all_y <- c(all_y, y_test[idx, 2])
    }

    graphics::plot(NA, xlim = range(all_x), ylim = range(all_y),
                   xlab = x$response_names[1], ylab = x$response_names[2],
                   main = sprintf("Localized Ellipsoids (alpha = %.2f)", x$alpha), ...)

    cols <- c("#2B83BA", "#D7191C", "#ABDDA4", "#FDAE61", "#984EA3", "#A65628")
    theta <- seq(0, 2 * pi, length.out = 100)
    circle <- rbind(cos(theta), sin(theta))

    for (ii in seq_along(idx)) {
      j <- idx[ii]
      mu <- x$pred[j, 1:2]
      Sigma <- x$cov_matrices[[j]][1:2, 1:2]
      q_val <- x$quantiles[j]

      eig <- eigen(Sigma)
      if (all(eig$values > 0)) {
        A <- eig$vectors %*% diag(sqrt(eig$values))
        ell <- t(mu + q_val * (A %*% circle))
        cl <- cols[(ii - 1) %% length(cols) + 1]
        graphics::polygon(ell[, 1], ell[, 2], border = cl, col = grDevices::adjustcolor(cl, alpha.f = 0.2), lwd = 2)
        graphics::points(mu[1], mu[2], pch = 4, col = cl, lwd = 2)
        graphics::text(mu[1], mu[2], labels = paste0("s_", j), pos = 3, cex = 0.8, col = cl)

        if (!is.null(y_test)) {
          y_obs <- y_test[j, 1:2]
          diff_y <- y_obs - mu
          dist_mah <- sqrt(as.numeric(crossprod(diff_y, solve(Sigma, diff_y))))
          is_cov <- (dist_mah <= q_val)
          graphics::points(y_obs[1], y_obs[2], pch = if (is_cov) 19 else 8,
                           col = if (is_cov) "darkgreen" else "red", cex = 1.2)
        }
      }
    }

    # Panel 2: Ellipsoid Volume across Space
    graphics::plot(x$s_test[, 1], x$ellipsoid_volumes,
                   xlab = "Spatial Coord 1", ylab = "Ellipsoid Volume",
                   main = "Spatial Uncertainty (Ellipsoid Volume)",
                   pch = 19, col = "#2B83BA", ...)
    graphics::lines(stats::lowess(x$s_test[, 1], x$ellipsoid_volumes), col = "red", lwd = 2)

  } else {
    mean_w <- colMeans(x$upper - x$lower)
    graphics::barplot(mean_w, names.arg = x$response_names, col = "steelblue",
                      ylab = "Mean Projected Width",
                      main = sprintf("Multivariate Conformal Widths (d = %d)", d), ...)
  }
  invisible(x)
}

#' Coverage Diagnostics for Multivariate Conformal Predictions
#'
#' Evaluates joint ellipsoidal containment and marginal coordinate coverage.
#'
#' @param object An object of class \code{"spconform_multivariate"}.
#' @param y_test Numeric matrix (\eqn{m \times d}) of true observed responses.
#'
#' @return A list with joint ellipsoidal coverage and marginal coverages.
#' @export
multivariate_coverage_report <- function(object, y_test) {
  if (!inherits(object, "spconform_multivariate")) {
    stop("object must be of class 'spconform_multivariate'.")
  }
  y_test <- as.matrix(y_test)
  m <- nrow(object$pred)
  d <- object$dimension

  if (nrow(y_test) != m || ncol(y_test) != d) {
    stop("y_test dimensions must match object predictions.")
  }

  joint_covered <- logical(m)
  marginal_covered <- matrix(FALSE, nrow = m, ncol = d, dimnames = list(NULL, object$response_names))

  for (j in seq_len(m)) {
    mu_j <- object$pred[j, ]
    Sigma_j <- object$cov_matrices[[j]]
    q_j <- object$quantiles[j]

    diff_j <- y_test[j, ] - mu_j
    dist_mah <- sqrt(as.numeric(crossprod(diff_j, solve(Sigma_j, diff_j))))
    joint_covered[j] <- (dist_mah <= q_j)

    for (k in seq_len(d)) {
      marginal_covered[j, k] <- (y_test[j, k] >= object$lower[j, k]) && (y_test[j, k] <= object$upper[j, k])
    }
  }

  res <- list(
    joint_coverage     = mean(joint_covered),
    marginal_coverage  = colMeans(marginal_covered),
    nominal_confidence = 1 - object$alpha,
    n_test             = m,
    dimension          = d,
    joint_covered_vec  = joint_covered
  )
  class(res) <- "spconform_multicov_report"
  res
}

#' Print Method for spconform_multicov_report
#' @param x An object of class \code{"spconform_multicov_report"}.
#' @param ... Additional arguments.
#' @export
print.spconform_multicov_report <- function(x, ...) {
  cat("=== Multivariate Conformal Coverage Evaluation ===\n")
  cat(sprintf("  Nominal Confidence Level : %.1f%%\n", 100 * x$nominal_confidence))
  cat(sprintf("  Empirical Joint Coverage : %.2f%% (Ellipsoidal Containment)\n", 100 * x$joint_coverage))
  cat(sprintf("  Sample Size (Test)       : %d locations\n", x$n_test))
  cat("\nMarginal Coordinate Coverage:\n")
  for (k in seq_along(x$marginal_coverage)) {
    nm <- names(x$marginal_coverage)[k]
    cat(sprintf("  Variable '%s': %.2f%%\n", nm, 100 * x$marginal_coverage[k]))
  }
  invisible(x)
}
