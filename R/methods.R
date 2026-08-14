#' @export
print.spconform <- function(x, ...) {
  cat("<spconform> ", x$type, " conformal prediction\n", sep = "")
  cat(sprintf("Target coverage: %.1f%%\n", 100 * (1 - x$alpha)))
  n <- length(x$pred)
  cat(sprintf("Number of prediction points: %d\n", n))
  k <- min(n, 6)
  df <- data.frame(
    pred = round(x$pred[seq_len(k)], 3),
    lower = round(x$lower[seq_len(k)], 3),
    upper = round(x$upper[seq_len(k)], 3)
  )
  print(df)
  if (n > k) cat(sprintf("... (%d more)\n", n - k))
  invisible(x)
}

#' @export
summary.spconform <- function(object, ...) {
  width <- object$upper - object$lower
  cat("spconform summary\n")
  cat("------------------\n")
  cat("Type:              ", object$type, "\n")
  cat("Target coverage:   ", sprintf("%.1f%%", 100 * (1 - object$alpha)), "\n")
  cat("Mean interval width:", round(mean(width), 4), "\n")
  cat("Median interval width:", round(stats::median(width), 4), "\n")
  invisible(object)
}


#' @importFrom graphics points segments legend
#' @export
plot.spconform <- function(x, y_true = NULL, ...) {
  n <- length(x$pred)
  ord <- seq_len(n)
  ylim <- range(c(x$lower, x$upper, x$pred, y_true), na.rm = TRUE)
  plot(ord, x$pred, ylim = ylim, pch = 19, xlab = "Index", ylab = "Predicted value",
       main = paste0("spconform (", x$type, ") - ",
                      round(100 * (1 - x$alpha)), "% intervals"), ...)
  segments(ord, x$lower, ord, x$upper, col = "grey50")
  if (!is.null(y_true)) {
    points(ord, y_true, col = "red", pch = 4)
    legend("topright", legend = c("prediction", "truth"),
           pch = c(19, 4), col = c("black", "red"), bty = "n")
  }
  invisible(x)
}

#' Empirical coverage and average interval width for an spconform object
#'
#' @param object An object of class \code{"spconform"}.
#' @param y_true Numeric vector of true observed values at the prediction
#'   locations, same length/order as \code{object$pred}.
#'
#' @return A named list with \code{coverage} (proportion of \code{y_true}
#'   falling within \code{[lower, upper]}) and \code{mean_width} (average
#'   interval width).
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' adj <- matrix(0, n, n)
#' for (i in 1:(n - 1)) { adj[i, i + 1] <- 1; adj[i + 1, i] <- 1 }
#' y <- cumsum(rnorm(n))
#' out <- scp_areal(y, adjacency = adj, alpha = 0.2)
#' coverage_report(out, y)
#'
#' @export
coverage_report <- function(object, y_true) {
  if (!inherits(object, "spconform")) stop("'object' must be of class 'spconform'.")
  covered <- (y_true >= object$lower) & (y_true <= object$upper)
  list(
    coverage = mean(covered),
    mean_width = mean(object$upper - object$lower)
  )
}
