#' Locally Weighted Split Conformal Prediction for Geostatistical (Point-Referenced) Data
#'
#' Constructs distribution-free prediction intervals at new spatial (or
#' spatio-temporal) locations by combining split conformal prediction with
#' spatial-distance kernel weights, relaxing the exchangeability assumption
#' that standard conformal prediction relies on.
#'
#' @param s_train Numeric matrix of training coordinates, one row per
#'   observation (e.g. longitude/latitude, or projected easting/northing).
#' @param y_train Numeric vector of training responses, same length as
#'   \code{nrow(s_train)}.
#' @param s0 Numeric matrix of coordinates at which prediction intervals are
#'   desired, one row per target location.
#' @param pred_fun A function with signature \code{function(s_train, y_train,
#'   s_new)} that returns a numeric vector of point predictions for the
#'   locations in \code{s_new}, fitted using \code{s_train}/\code{y_train}.
#'   This lets the user plug in any spatial predictor (kriging, random
#'   forest, GAM, etc.); \code{spconform} handles only the conformal
#'   calibration layer.
#' @param alpha Miscoverage level; prediction intervals target
#'   \eqn{1 - \alpha} coverage. Default 0.1 (90 percent intervals).
#' @param split Proportion of the training data used for model fitting; the
#'   remainder is used as the calibration set for conformal scores. Default
#'   0.5.
#' @param bandwidth Optional spatial kernel bandwidth passed to
#'   \code{\link{spatial_kernel_weights}}; if \code{NULL}, chosen
#'   automatically per target location.
#' @param t_train,t0,temporal_bandwidth Optional time indices (and temporal
#'   bandwidth) for spatio-temporal weighting; see
#'   \code{\link{spatial_kernel_weights}}.
#' @param seed Optional integer seed for the train/calibration split, for
#'   reproducibility.
#'
#' @return An object of class \code{"spconform"}, a list with components:
#'   \describe{
#'     \item{pred}{Point predictions at \code{s0}.}
#'     \item{lower, upper}{Lower/upper prediction interval bounds at
#'       \code{s0}.}
#'     \item{alpha}{The requested miscoverage level.}
#'     \item{s0}{The target coordinates.}
#'     \item{call}{The matched call.}
#'   }
#'
#' @details
#' This implements a localized split-conformal procedure in the spirit of
#' Mao, Martin and Reich (2020, "Valid model-free spatial prediction") and
#' subsequent spatio-temporal extensions: nonconformity scores from a
#' held-out calibration set are combined into a weighted empirical quantile
#' where weights decay with distance (in space, and optionally time) from
#' the target location, so nearby, more relevant observations dominate the
#' calibration of the interval while validity is retained under a local
#' exchangeability argument.
#'
#' @references
#' Mao, H., Martin, R., and Reich, B. J. (2020). Valid model-free spatial
#' prediction. \emph{arXiv:2006.15640}.
#'
#' @examples
#' set.seed(42)
#' n <- 200
#' s <- matrix(runif(2 * n), ncol = 2)
#' f <- function(s) sin(4 * s[, 1]) + cos(3 * s[, 2])
#' y <- f(s) + rnorm(n, sd = 0.3)
#'
#' idx <- sample(n, 150)
#' s_train <- s[idx, ]; y_train <- y[idx]
#' s0 <- s[-idx, ][1:10, ]
#'
#' pred_fun <- function(s_train, y_train, s_new) {
#'   fit <- lm(y_train ~ s_train[, 1] + s_train[, 2] +
#'               I(s_train[, 1]^2) + I(s_train[, 2]^2))
#'   nd <- data.frame(s_train = I(s_new))
#'   as.numeric(cbind(1, s_new[, 1], s_new[, 2], s_new[, 1]^2, s_new[, 2]^2) %*%
#'                coef(fit))
#' }
#'
#' out <- scp_geostatistical(s_train, y_train, s0, pred_fun, alpha = 0.1, seed = 1)
#' print(out)
#'
#' @export
scp_geostatistical <- function(s_train, y_train, s0, pred_fun,
                                alpha = 0.1, split = 0.5,
                                bandwidth = NULL,
                                t_train = NULL, t0 = NULL,
                                temporal_bandwidth = NULL,
                                seed = NULL) {
  cl <- match.call()

  s_train <- as.matrix(s_train)
  s0 <- as.matrix(s0)
  n <- nrow(s_train)

  if (length(y_train) != n) stop("'y_train' must match nrow(s_train).")
  if (!is.function(pred_fun)) stop("'pred_fun' must be a function.")
  if (alpha <= 0 || alpha >= 1) stop("'alpha' must be in (0, 1).")
  if (split <= 0 || split >= 1) stop("'split' must be in (0, 1).")

  if (!is.null(t_train) && length(t_train) != n) {
    stop("'t_train' must match nrow(s_train).")
  }
  if (!is.null(t0) && nrow(s0) != length(t0)) {
    stop("'t0' must have length equal to nrow(s0).")
  }

  if (!is.null(seed)) set.seed(seed)

  n_fit <- floor(split * n)
  fit_idx <- sample(seq_len(n), n_fit)
  cal_idx <- setdiff(seq_len(n), fit_idx)
  if (length(cal_idx) < 5) {
    stop("Too few calibration points; increase n or decrease 'split'.")
  }

  s_fit <- s_train[fit_idx, , drop = FALSE]
  y_fit <- y_train[fit_idx]
  s_cal <- s_train[cal_idx, , drop = FALSE]
  y_cal <- y_train[cal_idx]
  t_cal <- if (!is.null(t_train)) t_train[cal_idx] else NULL

  cal_pred <- pred_fun(s_fit, y_fit, s_cal)
  resid_abs <- abs(y_cal - cal_pred)

  pred0 <- pred_fun(s_fit, y_fit, s0)

  m <- nrow(s0)
  lower <- upper <- numeric(m)

  for (j in seq_len(m)) {
    t0j <- if (!is.null(t0)) t0[j] else NULL
    w <- spatial_kernel_weights(
      s0 = s0[j, ], s = s_cal, bandwidth = bandwidth,
      t0 = t0j, t = t_cal, temporal_bandwidth = temporal_bandwidth
    )

    ord <- order(resid_abs)
    r_sorted <- resid_abs[ord]
    w_sorted <- w[ord]
    w_sorted <- w_sorted + 1e-12
    cw <- cumsum(w_sorted) / sum(w_sorted)

    q_level <- min(1, (1 - alpha) * (length(r_sorted) + 1) / length(r_sorted))
    q_idx <- which(cw >= q_level)[1]
    if (is.na(q_idx)) q_idx <- length(r_sorted)
    qhat <- r_sorted[q_idx]

    lower[j] <- pred0[j] - qhat
    upper[j] <- pred0[j] + qhat
  }

  structure(
    list(
      pred = pred0,
      lower = lower,
      upper = upper,
      alpha = alpha,
      s0 = s0,
      type = "geostatistical",
      call = cl
    ),
    class = "spconform"
  )
}
