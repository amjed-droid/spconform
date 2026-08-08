#' Compute Gaussian-kernel spatial weights between a target location and a
#' set of reference locations
#'
#' @param s0 A numeric vector (or 1-row matrix) giving the coordinates of the
#'   target (prediction) location.
#' @param s A numeric matrix of coordinates (one row per observation) for the
#'   reference (calibration) set.
#' @param bandwidth Positive numeric bandwidth of the Gaussian kernel. If
#'   \code{NULL} (default), the bandwidth is chosen automatically as the
#'   median pairwise distance among \code{s} (Silverman-type rule).
#' @param t0 Optional numeric scalar: time index of the target observation,
#'   for spatio-temporal weighting. If supplied, \code{t} must be supplied too.
#' @param t Optional numeric vector: time indices for the reference set,
#'   same length as \code{nrow(s)}.
#' @param temporal_bandwidth Positive numeric bandwidth for the temporal
#'   kernel; only used if \code{t0}/\code{t} are supplied. Defaults to the
#'   spatial bandwidth logic applied to \code{t}.
#'
#' @return A numeric vector of weights (not necessarily summing to one),
#'   one per row of \code{s}, giving higher weight to reference points that
#'   are spatially (and, if requested, temporally) closer to \code{s0}.
#'
#' @details Weights are computed as
#'   \deqn{w_i = \exp(-\|s_0 - s_i\|^2 / (2 h^2))}
#'   for the spatial-only case, and multiplied by an analogous temporal
#'   kernel term when \code{t0}/\code{t} are provided. This is the core
#'   device used to relax the exchangeability assumption required by
#'   standard conformal prediction: observations located near the point to
#'   be predicted contribute more to the calibration of the prediction
#'   interval than distant ones (Mao, Martin and Reich, 2020).
#'
#' @examples
#' set.seed(1)
#' s <- matrix(runif(20), ncol = 2)
#' s0 <- c(0.5, 0.5)
#' w <- spatial_kernel_weights(s0, s)
#' round(w, 3)
#'
#' @export
spatial_kernel_weights <- function(s0, s, bandwidth = NULL,
                                    t0 = NULL, t = NULL,
                                    temporal_bandwidth = NULL) {
  s <- as.matrix(s)
  s0 <- as.numeric(s0)
  if (length(s0) != ncol(s)) {
    stop("'s0' must have the same number of coordinates as columns in 's'.")
  }

  d <- sqrt(rowSums(sweep(s, 2, s0, "-")^2))

  if (is.null(bandwidth)) {
    pd <- stats::dist(s)
    bandwidth <- stats::median(pd)
    if (!is.finite(bandwidth) || bandwidth <= 0) bandwidth <- 1
  }

  w <- exp(-(d^2) / (2 * bandwidth^2))

  if (!is.null(t0) || !is.null(t)) {
    if (is.null(t0) || is.null(t)) {
      stop("Both 't0' and 't' must be supplied for spatio-temporal weighting.")
    }
    if (length(t) != nrow(s)) {
      stop("'t' must have length equal to nrow(s).")
    }
    if (is.null(temporal_bandwidth)) {
      td <- abs(outer(t, t, "-"))
      temporal_bandwidth <- stats::median(td[upper.tri(td)])
      if (!is.finite(temporal_bandwidth) || temporal_bandwidth <= 0) {
        temporal_bandwidth <- 1
      }
    }
    wt <- exp(-((t - t0)^2) / (2 * temporal_bandwidth^2))
    w <- w * wt
  }

  w
}

#' Compute neighbourhood-based weights for areal (lattice) data
#'
#' @param i0 Integer index of the target (unobserved / held-out) areal unit.
#' @param adjacency A square 0/1 (or weighted) adjacency matrix describing
#'   the neighbourhood structure of the areal units (e.g. a spatial contiguity
#'   matrix). Row/column \code{i0} corresponds to the target unit.
#' @param decay Numeric decay rate applied to graph distance (number of hops)
#'   from \code{i0}; larger values down-weight distant neighbours more
#'   aggressively. Defaults to 1.
#'
#' @return A numeric vector of length \code{nrow(adjacency)} with weights
#'   based on graph distance from \code{i0} (self-weight is 0, i.e. the
#'   target unit is excluded from its own calibration set).
#'
#' @export
areal_neighbor_weights <- function(i0, adjacency, decay = 1) {
  adjacency <- as.matrix(adjacency)
  n <- nrow(adjacency)
  if (n != ncol(adjacency)) stop("'adjacency' must be a square matrix.")
  if (i0 < 1 || i0 > n) stop("'i0' out of range.")

  adj_bin <- (adjacency != 0) * 1
  diag(adj_bin) <- 0

  dist_hops <- rep(Inf, n)
  dist_hops[i0] <- 0
  frontier <- i0
  hop <- 0
  visited <- rep(FALSE, n)
  visited[i0] <- TRUE
  while (length(frontier) > 0) {
    hop <- hop + 1
    nbrs <- which(colSums(adj_bin[frontier, , drop = FALSE] != 0) > 0)
    nbrs <- nbrs[!visited[nbrs]]
    if (length(nbrs) == 0) break
    dist_hops[nbrs] <- hop
    visited[nbrs] <- TRUE
    frontier <- nbrs
  }

  w <- exp(-decay * dist_hops)
  w[i0] <- 0
  w[!is.finite(w)] <- 0
  w
}
