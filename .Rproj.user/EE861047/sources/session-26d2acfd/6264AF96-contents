#' Neighbourhood-Weighted Conformal Prediction for Areal (Lattice) Data
#'
#' Constructs distribution-free prediction intervals for areal units (e.g.
#' counties, census tracts, grid cells) using a leave-one-unit-out
#' conformal procedure in which nonconformity scores are weighted by graph
#' (neighbourhood) distance to the held-out unit, via its adjacency
#' structure.
#'
#' @param y Numeric vector of observed responses, one per areal unit.
#' @param X Optional numeric design matrix of covariates (one row per areal
#'   unit); passed to \code{pred_fun} if supplied. May be \code{NULL} for
#'   purely spatial (neighbourhood-based) prediction.
#' @param adjacency Square adjacency (contiguity) matrix describing the
#'   neighbourhood structure among the \code{length(y)} areal units.
#'   Non-zero entries are treated as neighbours (coerced to binary).
#' @param pred_fun A function with signature \code{function(y_train, X_train,
#'   idx_train, idx_target, adjacency)} returning a single numeric point
#'   prediction for the areal unit indexed by \code{idx_target}, fitted
#'   using the training units \code{idx_train}. If \code{NULL} (default), a
#'   simple neighbourhood-mean predictor is used (the mean of \code{y_train}
#'   over the target's graph neighbours, falling back to the global training
#'   mean if the target has no observed neighbours).
#' @param alpha Miscoverage level; intervals target \eqn{1-\alpha} coverage.
#'   Default 0.1.
#' @param decay Decay rate for neighbourhood weights; see
#'   \code{\link{areal_neighbor_weights}}.
#'
#' @return An object of class \code{"spconform"} (areal variant) with the
#'   same structure as \code{\link{scp_geostatistical}}, using unit indices
#'   in place of coordinates.
#'
#' @details Uses a full leave-one-out (jackknife-style) conformal scheme:
#'   for each areal unit \eqn{i}, a model is fit on all other units and used
#'   to predict unit \eqn{i}; the resulting nonconformity scores across all
#'   units are combined into a graph-distance-weighted quantile specific to
#'   each target unit, so that the calibration set is dominated by
#'   spatially/graph-proximate units rather than treating all units as
#'   exchangeable.
#'
#'   If the effective neighbourhood weight for a target unit is
#'   vanishingly small (e.g., an isolated unit with no neighbours), the
#'   procedure falls back to an unweighted quantile over all other units
#'   and issues a warning.
#'
#' @examples
#' set.seed(1)
#' n <- 25
#' adj <- matrix(0, n, n)
#' for (i in 1:(n - 1)) { adj[i, i + 1] <- 1; adj[i + 1, i] <- 1 }
#' y <- cumsum(rnorm(n)) + rnorm(n, sd = 0.2)
#'
#' out <- scp_areal(y, adjacency = adj, alpha = 0.2)
#' print(out)
#'
#' @export
scp_areal <- function(y, X = NULL, adjacency, pred_fun = NULL,
                      alpha = 0.1, decay = 1) {
  cl <- match.call()
  
  n <- length(y)
  adjacency <- as.matrix(adjacency)
  if (nrow(adjacency) != n || ncol(adjacency) != n) {
    stop("'adjacency' must be an n x n matrix matching length(y).")
  }
  if (alpha <= 0 || alpha >= 1) stop("'alpha' must be in (0, 1).")
  
  if (!is.null(X)) {
    X <- as.matrix(X)
    if (nrow(X) != n) {
      stop("'X' must have same number of rows as length(y).")
    }
  }
  
  if (any(adjacency < 0, na.rm = TRUE)) {
    stop("'adjacency' must not contain negative values.")
  }
  if (any(adjacency > 1, na.rm = TRUE)) {
    warning("'adjacency' contains values > 1; coercing to binary (non-zero -> 1).")
    adjacency <- (adjacency != 0) * 1
  }
  diag(adjacency) <- 0
  
  if (is.null(pred_fun)) {
    pred_fun <- function(y_train, X_train, idx_train, idx_target, adjacency) {
      adj_bin <- (adjacency != 0) * 1
      nbrs <- which(adj_bin[idx_target, ] != 0)
      nbrs <- intersect(nbrs, idx_train)
      if (length(nbrs) == 0) {
        return(mean(y_train))
      }
      # Convert original indices to relative indices within y_train
      rel_idx <- match(nbrs, idx_train)
      mean(y_train[rel_idx])
    }
  }
  
  loo_pred <- numeric(n)
  for (i in seq_len(n)) {
    idx_train <- setdiff(seq_len(n), i)
    X_train <- if (!is.null(X)) X[idx_train, , drop = FALSE] else NULL
    loo_pred[i] <- pred_fun(y[idx_train], X_train, idx_train, i, adjacency)
  }
  resid_abs <- abs(y - loo_pred)
  
  pred <- loo_pred
  lower <- upper <- numeric(n)
  
  for (i in seq_len(n)) {
    w <- areal_neighbor_weights(i, adjacency, decay = decay)
    others <- setdiff(seq_len(n), i)
    
    r_i <- resid_abs[others]
    # التأكد من خلوه من القيم المفقودة
    if (anyNA(r_i)) {
      warning(sprintf("Nonconformity scores contain NA for unit %d; using unweighted quantile.", i))
      w_i <- rep(1, length(others))
    } else {
      w_i <- w[others] + 1e-12
      # إذا كانت الأوزان الفعلية كلها صفراً تقريباً، استخدم أوزاناً متساوية
      if (sum(w_i) < 1e-10 || max(w_i) < 1e-12) {
        warning(sprintf("Unit %d has negligible neighbourhood weights; falling back to unweighted quantile.", i))
        w_i <- rep(1, length(others))
      }
    }
    
    ord <- order(r_i)
    r_sorted <- r_i[ord]
    w_sorted <- w_i[ord]
    cw <- cumsum(w_sorted) / sum(w_sorted)
    
    q_level <- min(1, (1 - alpha) * (length(r_sorted) + 1) / length(r_sorted))
    q_idx <- which(cw >= q_level)[1]
    if (is.na(q_idx)) {
      # إذا لم يتم العثور على كمّاء، استخدم أكبر قيمة (حالة متطرفة)
      q_idx <- length(r_sorted)
    }
    qhat <- r_sorted[q_idx]
    
    lower[i] <- pred[i] - qhat
    upper[i] <- pred[i] + qhat
  }
  
  structure(
    list(
      pred = pred,
      lower = lower,
      upper = upper,
      alpha = alpha,
      s0 = seq_len(n),
      type = "areal",
      call = cl
    ),
    class = "spconform"
  )
}