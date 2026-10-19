#' Spatially Localized Conformal Prediction for Classification
#'
#' Constructs distribution-free, model-agnostic prediction sets for binary and
#' multiclass spatial data under spatial autocorrelation using localized conformal
#' calibration (Adaptive Prediction Sets with spatial kernel weighting).
#'
#' @param s_train Numeric matrix or data frame of training coordinates (N_train x 2).
#' @param y_train Factor or character vector of training class labels (length N_train).
#' @param s0 Numeric matrix or data frame of target coordinates for prediction (N_test x 2).
#' @param pred_fun Function with signature \code{function(s_train, y_train, s_new)} that
#'   fits a classifier and returns an N_new x K matrix of predicted probabilities whose
#'   rows sum to 1 and column names match the class labels.
#' @param alpha Miscoverage level; prediction sets guarantee at least \eqn{1 - \alpha}
#'   spatial coverage. Default: 0.10 (90 percent confidence).
#' @param split Proportion of training data used for model fitting; the remainder
#'   is used as the calibration set for conformal score evaluation. Default: 0.50.
#' @param bandwidth Optional numeric spatial kernel bandwidth. Defaults to the median
#'   pairwise distance of calibration points.
#' @param method Nonconformity scoring method: \code{"APS"} (Adaptive Prediction Sets,
#'   cumulative probability, default) or \code{"LAC"} (Least Likely / Top-class score).
#'
#' @return An object of class \code{"spconform_classification"} containing:
#' \describe{
#'   \item{prediction_sets}{A list of length N_test containing character vectors of
#'     admissible class labels for each location.}
#'   \item{set_sizes}{Integer vector giving the number of classes in each prediction set.}
#'   \item{prob_pred}{Matrix of predicted class probabilities at test locations.}
#'   \item{classes}{Character vector of unique class labels.}
#'   \item{alpha}{Nominal miscoverage rate.}
#'   \item{bandwidth}{Spatial kernel bandwidth used.}
#'   \item{s0}{Target spatial coordinates matrix.}
#'   \item{method}{Score formulation used.}
#' }
#' @export
#'
#' @examples
#' set.seed(123)
#' n_tr <- 100
#' s_tr <- matrix(runif(n_tr * 2), n_tr, 2)
#' y_tr <- factor(ifelse(s_tr[, 1] + s_tr[, 2] > 1.0, "Class_B", "Class_A"))
#' s_te <- matrix(runif(30 * 2), 30, 2)
#' 
#' # Define simple spatial logistic regression predictor
#' pred_fun <- function(s_train, y_train, s_new) {
#'   df_tr <- data.frame(y = y_train, x = s_train[, 1], y_coord = s_train[, 2])
#'   df_nw <- data.frame(x = s_new[, 1], y_coord = s_new[, 2])
#'   fit <- glm(y ~ x + y_coord, data = df_tr, family = binomial)
#'   p_b <- predict(fit, newdata = df_nw, type = "response")
#'   cbind(Class_A = 1 - p_b, Class_B = p_b)
#' }
#' 
#' fit_class <- scp_classification(s_tr, y_tr, s_te, pred_fun, alpha = 0.10)
#' print(fit_class)
scp_classification <- function(s_train, y_train, s0, pred_fun, 
                               alpha = 0.10, split = 0.50, 
                               bandwidth = NULL,
                               method = c("APS", "LAC")) {
  method <- match.arg(method)
  cl <- match.call()
  
  if (alpha <= 0 || alpha >= 1) stop("'alpha' must be strictly between 0 and 1.")
  if (split <= 0 || split >= 1) stop("'split' must be strictly between 0 and 1.")
  if (!is.function(pred_fun)) stop("'pred_fun' must be a valid prediction function.")
  
  s_train <- as.matrix(s_train)
  s0      <- as.matrix(s0)
  n_total <- nrow(s_train)
  
  if (n_total < 10) stop("Training data requires at least n = 10 observations.")
  if (length(y_train) != n_total) stop("length(y_train) must equal nrow(s_train).")
  
  # Ensure y_train is a factor
  y_train <- as.factor(y_train)
  classes <- levels(y_train)
  K <- length(classes)
  if (K < 2) stop("Classification requires at least 2 distinct classes.")
  
  # Split training data into proper training and calibration subsets
  n_fit <- max(2, round(n_total * split))
  idx_fit <- sample.int(n_total, size = n_fit)
  idx_cal <- setdiff(seq_len(n_total), idx_fit)
  n_cal <- length(idx_cal)
  
  s_fit <- s_train[idx_fit, , drop = FALSE]
  y_fit <- y_train[idx_fit]
  s_cal <- s_train[idx_cal, , drop = FALSE]
  y_cal <- y_train[idx_cal]
  
  # Fit base classifier and predict probabilities on calibration and target sets
  prob_cal <- as.matrix(pred_fun(s_fit, y_fit, s_cal))
  prob_0   <- as.matrix(pred_fun(s_fit, y_fit, s0))
  
  # Ensure column names match class labels
  if (is.null(colnames(prob_cal))) colnames(prob_cal) <- classes
  if (is.null(colnames(prob_0)))   colnames(prob_0)   <- classes
  
  # Normalize rows to sum strictly to 1
  prob_cal <- prob_cal / rowSums(prob_cal)
  prob_0   <- prob_0 / rowSums(prob_0)
  
  # Evaluate nonconformity scores on calibration set
  # E_i for each calibration observation
  cal_scores <- numeric(n_cal)
  for (i in seq_len(n_cal)) {
    true_cls <- as.character(y_cal[i])
    p_vec <- prob_cal[i, ]
    
    if (method == "LAC") {
      # Least Likely Score: E_i = 1 - p(Y_i)
      cal_scores[i] <- 1 - p_vec[true_cls]
    } else {
      # Adaptive Prediction Sets (APS): cumulative probability of sorted classes down to true label
      ord <- order(p_vec, decreasing = TRUE)
      sorted_cls <- names(p_vec)[ord]
      sorted_p   <- p_vec[ord]
      k_star <- which(sorted_cls == true_cls)[1]
      cal_scores[i] <- sum(sorted_p[seq_len(k_star)])
    }
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
  
  # Spatially localized prediction sets for each target location
  n_test <- nrow(s0)
  prediction_sets <- vector("list", n_test)
  set_sizes <- integer(n_test)
  
  for (j in seq_len(n_test)) {
    # Spatial kernel weights: w_i(s0_j)
    w_vec <- exp(- (dist_cal_0[, j]^2) / (2 * bandwidth^2)) + 1e-12
    w_norm <- w_vec / sum(w_vec)
    
    # Sort calibration scores and compute spatially weighted quantile
    ord_sc <- order(cal_scores)
    scores_sorted <- cal_scores[ord_sc]
    w_sorted <- w_norm[ord_sc]
    cum_w <- cumsum(w_sorted)
    
    q_level <- min(1.0, (1 - alpha) * (1 + 1 / n_cal))
    idx_q <- which(cum_w >= q_level)[1]
    if (is.na(idx_q)) idx_q <- length(scores_sorted)
    q_val <- scores_sorted[idx_q]
    
    # Form prediction set at s0_j
    p_test_j <- prob_0[j, ]
    ord_test <- order(p_test_j, decreasing = TRUE)
    cls_sorted <- names(p_test_j)[ord_test]
    probs_sorted <- p_test_j[ord_test]
    
    if (method == "LAC") {
      # LAC: include class k if 1 - p_k <= q_val (i.e. p_k >= 1 - q_val)
      threshold <- 1 - q_val
      admissible <- names(p_test_j)[p_test_j >= threshold]
      if (length(admissible) == 0) {
        # Fallback to top-1 class if empty
        admissible <- cls_sorted[1]
      }
    } else {
      # APS (Adaptive Prediction Sets): include classes while the cumulative
      # probability of previous classes does not exceed the calibrated threshold q_val
      cum_prev <- c(0, cumsum(probs_sorted[-length(probs_sorted)]))
      admissible <- cls_sorted[cum_prev <= q_val]
      if (length(admissible) == 0) admissible <- cls_sorted[1]
    }
    
    prediction_sets[[j]] <- admissible
    set_sizes[j] <- length(admissible)
  }
  
  structure(
    list(
      prediction_sets = prediction_sets,
      set_sizes = set_sizes,
      prob_pred = prob_0,
      classes = classes,
      alpha = alpha,
      bandwidth = round(bandwidth, 4),
      s0 = s0,
      method = method,
      call = cl
    ),
    class = "spconform_classification"
  )
}

#' @export
print.spconform_classification <- function(x, ...) {
  cat("-----------------------------------------------------------\n")
  cat("  <spconform> Spatial Conformal Classification\n")
  cat("-----------------------------------------------------------\n")
  cat(sprintf("Target Nominal Coverage:      %.1f%%\n", 100 * (1 - x$alpha)))
  cat("Scoring Method:               ", x$method, "\n", sep = "")
  cat(sprintf("Number of Prediction Points:  %d\n", length(x$set_sizes)))
  cat("Classes Evaluated:            ", paste(x$classes, collapse = ", "), "\n", sep = "")
  cat(sprintf("Spatial Bandwidth (h):        %.4f\n", x$bandwidth))
  cat("-----------------------------------------------------------\n")
  cat(sprintf("Average Prediction Set Size:  %.2f classes\n", mean(x$set_sizes)))
  cat(sprintf("Singleton Rate (Definitive):  %.1f%%\n", 100 * mean(x$set_sizes == 1)))
  cat(sprintf("Ambiguous Sets (> 1 class):   %.1f%%\n", 100 * mean(x$set_sizes > 1)))
  cat("-----------------------------------------------------------\n")
  invisible(x)
}

#' @export
summary.spconform_classification <- function(object, ...) {
  cat("===========================================================\n")
  cat("  spconform: Spatial Classification Summary\n")
  cat("===========================================================\n")
  print(object)
  cat("\nSet Size Distribution:\n")
  tbl <- table(object$set_sizes)
  pcts <- prop.table(tbl) * 100
  for (sz in names(tbl)) {
    cat(sprintf("  Size %s : %3d locations (%5.1f%%)\n", sz, tbl[sz], pcts[sz]))
  }
  cat("===========================================================\n")
  invisible(object)
}

#' @export
plot.spconform_classification <- function(x, ...) {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  
  graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  
  # Colors for set size: 1 = Blue/Green (Certain), 2 = Orange (Ambiguous), 3+ = Red (High Uncertainty)
  size_cols <- ifelse(x$set_sizes == 1, "#2B83BA",
               ifelse(x$set_sizes == 2, "#FDAE61", "#D7191C"))
  
  # Panel 1: Map of Prediction Set Sizes (Decisiveness Map)
  plot(x$s0[, 1], x$s0[, 2], col = size_cols, pch = 19, cex = 1.3,
       xlab = "Easting (X)", ylab = "Northing (Y)",
       main = "Prediction Set Size Map")
  graphics::grid(col = "gray85")
  graphics::points(x$s0[, 1], x$s0[, 2], col = "gray30", pch = 1, cex = 1.3)
  
  legend_labels <- c("Singleton (1 class)", "Ambiguous (2 classes)", "Uncertain (3+ classes)")
  legend_colors <- c("#2B83BA", "#FDAE61", "#D7191C")
  graphics::legend("topleft", legend = legend_labels, col = legend_colors,
                   pch = 19, pt.cex = 1.1, bty = "n", cex = 0.8)
  
  # Panel 2: Distribution of Set Sizes (Barplot)
  tbl <- table(factor(x$set_sizes, levels = 1:max(x$set_sizes)))
  bp <- graphics::barplot(tbl, col = legend_colors[1:length(tbl)], border = "gray30",
                          xlab = "Prediction Set Size", ylab = "Frequency",
                          main = "Set Size Distribution")
  graphics::grid(nx = NA, ny = NULL, col = "gray85")
  graphics::text(bp, tbl, labels = as.numeric(tbl), pos = 3, cex = 0.9)
}
