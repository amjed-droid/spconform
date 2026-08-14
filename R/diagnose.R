#' Comprehensive Diagnostic Report for spconform Objects
#'
#' Produces a multi-panel diagnostic report assessing marginal coverage,
#' conditional coverage by spatial strata, boundary effects, and the
#' distribution of nonconformity scores.
#'
#' @param object An object of class \code{"spconform"}, typically the output
#'   of \code{scp_geostatistical()} or \code{scp_areal()}.
#' @param y_true Numeric vector of true response values at the prediction
#'   locations. Must have the same length as \code{object$pred}.
#' @param s_test Optional numeric matrix of prediction coordinates
#'   (one row per location). Required for spatial diagnostics.
#' @param n_bins Integer; number of spatial bins for conditional coverage
#'   (default 4).
#' @param plot Logical; if \code{TRUE} (default), generates a multi-panel
#'   diagnostic plot.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return A list (invisibly) containing:
#' \describe{
#'   \item{marginal}{Marginal coverage and mean width.}
#'   \item{conditional}{Coverage and width by spatial bin.}
#'   \item{boundary}{Coverage by distance from convex hull boundary.}
#'   \item{scores}{Summary of nonconformity score distribution.}
#' }
#'
#' @details
#' The conditional coverage analysis partitions the prediction locations
#' into \code{n_bins} equal-area spatial quadrants and reports coverage
#' within each. The boundary analysis classifies points by their distance
#' to the convex hull of the training data (for geostatistical output).
#'
#' If empirical coverage exceeds the nominal target by more than 5
#' percentage points, a message is issued suggesting bandwidth
#' reduction for tighter intervals.
#'
#' @examples
#' \dontrun{
#' out <- scp_geostatistical(s_train, y_train, s_test, pred_fun, alpha = 0.1)
#' diag <- diagnose(out, y_test, s_test)
#' }
#'
#' @importFrom stats median sd quantile qqnorm qqline
#' @importFrom graphics abline barplot layout par text
#' @importFrom grDevices chull
#' @export
diagnose <- function(object, y_true, s_test = NULL, n_bins = 4, plot = TRUE, ...) {
  
  if (!inherits(object, "spconform")) {
    stop("'object' must be of class 'spconform'")
  }
  
  n <- length(object$pred)
  if (length(y_true) != n) {
    stop("length(y_true) must equal length(object$pred)")
  }
  
  ## ---- 1. Marginal coverage ----
  covered <- (y_true >= object$lower) & (y_true <= object$upper)
  widths  <- object$upper - object$lower
  
  marginal <- list(
    coverage    = mean(covered, na.rm = TRUE),
    mean_width  = mean(widths, na.rm = TRUE),
    median_width = median(widths, na.rm = TRUE),
    sd_width    = sd(widths, na.rm = TRUE),
    n           = n,
    n_covered   = sum(covered, na.rm = TRUE)
  )
  
  ## ---- Overcoverage warning ----
  nominal <- 1 - object$alpha
  if (marginal$coverage > nominal + 0.05) {
    message("Note: Empirical coverage (", round(marginal$coverage, 3),
            ") exceeds nominal (", nominal, ") by >5%.",
            " Consider reducing 'bandwidth' for tighter intervals.")
  } else if (marginal$coverage < nominal - 0.05) {
    warning("Empirical coverage (", round(marginal$coverage, 3),
            ") is >5% below nominal (", nominal, ").",
            " Consider increasing 'bandwidth' or checking model fit.")
  }
  
  ## ---- 2. Conditional coverage by spatial bins ----
  conditional <- NULL
  if (!is.null(s_test) && ncol(s_test) >= 2) {
    x_q <- cut(s_test[, 1], breaks = n_bins, include.lowest = TRUE, labels = FALSE)
    y_q <- cut(s_test[, 2], breaks = n_bins, include.lowest = TRUE, labels = FALSE)
    bin_id <- paste0("Q", x_q, "-", y_q)
    
    conditional <- data.frame(
      bin   = character(),
      n     = integer(),
      coverage = numeric(),
      mean_width = numeric(),
      stringsAsFactors = FALSE
    )
    
    for (b in sort(unique(bin_id))) {
      idx <- which(bin_id == b)
      conditional <- rbind(conditional, data.frame(
        bin        = b,
        n          = length(idx),
        coverage   = mean(covered[idx], na.rm = TRUE),
        mean_width = mean(widths[idx], na.rm = TRUE)
      ))
    }
  }
  
  ## ---- 3. Boundary effects ----
  boundary <- NULL
  if (!is.null(s_test) && ncol(s_test) >= 2 && object$type == "geostatistical") {
    hull <- chull(s_test[, 1], s_test[, 2])
    hull_pts <- s_test[hull, , drop = FALSE]
    
    dist_to_boundary <- apply(s_test, 1, function(pt) {
      min(sqrt((hull_pts[, 1] - pt[1])^2 + (hull_pts[, 2] - pt[2])^2))
    })
    
    boundary_thresh <- median(dist_to_boundary, na.rm = TRUE)
    is_near <- dist_to_boundary <= boundary_thresh
    
    boundary <- list(
      threshold   = boundary_thresh,
      near_boundary = list(
        n        = sum(is_near),
        coverage = mean(covered[is_near], na.rm = TRUE),
        mean_width = mean(widths[is_near], na.rm = TRUE)
      ),
      far_boundary = list(
        n        = sum(!is_near),
        coverage = mean(covered[!is_near], na.rm = TRUE),
        mean_width = mean(widths[!is_near], na.rm = TRUE)
      )
    )
  }
  
  ## ---- 4. Nonconformity score distribution ----
  scores <- pmax(object$upper - object$pred, object$pred - object$lower)
  
  scores_summary <- list(
    mean   = mean(scores, na.rm = TRUE),
    median = median(scores, na.rm = TRUE),
    sd     = sd(scores, na.rm = TRUE),
    q90    = quantile(scores, 0.90, na.rm = TRUE),
    q95    = quantile(scores, 0.95, na.rm = TRUE)
  )
  
  ## ---- 5. Plot ----
  if (plot) {
    n_panels <- 2
    if (!is.null(conditional)) n_panels <- n_panels + 1
    if (!is.null(boundary)) n_panels <- n_panels + 1
    
    layout_mat <- matrix(1:n_panels, ncol = 2, byrow = TRUE)
    if (n_panels %% 2 == 1) {
      layout_mat <- rbind(layout_mat, c(n_panels, n_panels))
    }
    
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    
    layout(layout_mat)
    par(mar = c(4, 4, 3, 1))
    
    ## Panel 1: Marginal coverage bar
    bp <- barplot(c(marginal$coverage, nominal),
                  names.arg = c("Empirical", "Nominal"),
                  col = c("steelblue", "gray80"),
                  ylab = "Coverage",
                  main = paste0("Marginal coverage (n=", marginal$n, ")"),
                  ylim = c(0, 1.15)) # زيادة بسيطة ليتسع النص بالأعلى
    abline(h = nominal, col = "red", lty = 2)
    
    ## وضع القيمة الرقمية بدقة فوق العمود الأول
    text(bp[1], min(marginal$coverage + 0.06, 1.06), 
         labels = round(marginal$coverage, 3), cex = 1.1, font = 2)
    ## وضع القيمة الاسمية فوق العمود الثاني
    text(bp[2], nominal + 0.06, 
         labels = round(nominal, 3), cex = 1.1, font = 2)
    
    ## Panel 2: Nonconformity score QQ-plot
    qqnorm(scores, main = "QQ-plot of nonconformity scores", 
           pch = 20, col = "steelblue")
    qqline(scores, col = "red", lty = 2)
    
    ## Panel 3: Conditional coverage by bin
    if (!is.null(conditional)) {
      barplot(conditional$coverage, names.arg = conditional$bin,
              col = "steelblue", ylab = "Coverage",
              main = "Conditional coverage by spatial bin",
              ylim = c(0, 1), las = 2, cex.names = 0.7)
      abline(h = nominal, col = "red", lty = 2)
    }
    
    ## Panel 4: Boundary effect
    if (!is.null(boundary)) {
      barplot(c(boundary$near_boundary$coverage, boundary$far_boundary$coverage),
              names.arg = c("Near boundary", "Far from boundary"),
              col = c("coral", "steelblue"), ylab = "Coverage",
              main = "Coverage vs boundary proximity",
              ylim = c(0, 1))
      abline(h = nominal, col = "red", lty = 2)
    }
  }
  
  ## ---- Return ----
  out <- list(
    marginal    = marginal,
    conditional = conditional,
    boundary    = boundary,
    scores      = scores_summary,
    alpha       = object$alpha
  )
  
  class(out) <- c("spconform_diagnose", "list")
  invisible(out)
}


#' @export
print.spconform_diagnose <- function(x, ...) {
  cat("=== spconform Diagnostic Report ===\n\n")
  
  cat("Marginal coverage:\n")
  cat("  Empirical:", round(x$marginal$coverage, 4), 
      "(nominal:", 1 - x$alpha, ")\n")
  cat("  Mean width:", round(x$marginal$mean_width, 4), "\n")
  cat("  n =", x$marginal$n, ", covered =", x$marginal$n_covered, "\n\n")
  
  if (!is.null(x$conditional)) {
    cat("Conditional coverage by spatial bin:\n")
    for (i in seq_len(nrow(x$conditional))) {
      cat("  ", x$conditional$bin[i], ": ", 
          round(x$conditional$coverage[i], 4),
          " (n=", x$conditional$n[i], ", width=", 
          round(x$conditional$mean_width[i], 3), ")\n", sep = "")
    }
    cat("\n")
  }
  
  if (!is.null(x$boundary)) {
    cat("Boundary effect:\n")
    cat("  Near boundary:  ", round(x$boundary$near_boundary$coverage, 4),
        " (n=", x$boundary$near_boundary$n, ")\n", sep = "")
    cat("  Far from boundary:", round(x$boundary$far_boundary$coverage, 4),
        " (n=", x$boundary$far_boundary$n, ")\n\n", sep = "")
  }
  
  cat("Nonconformity scores:\n")
  cat("  Mean:", round(x$scores$mean, 4), "\n")
  cat("  Median:", round(x$scores$median, 4), "\n")
  cat("  SD:", round(x$scores$sd, 4), "\n")
  cat("  90% quantile:", round(x$scores$q90, 4), "\n")
}