test_that("spatial_kernel_weights returns valid weights", {
  set.seed(1)
  s <- matrix(runif(20), ncol = 2)
  s0 <- c(0.5, 0.5)
  w <- spatial_kernel_weights(s0, s)
  expect_length(w, nrow(s))
  expect_true(all(w >= 0 & w <= 1))
})

test_that("scp_geostatistical achieves roughly nominal coverage", {
  set.seed(42)
  n <- 400
  s <- matrix(runif(2 * n), ncol = 2)
  f <- function(s) sin(4 * s[, 1]) + cos(3 * s[, 2])
  y <- f(s) + rnorm(n, sd = 0.3)

  idx <- sample(n, 300)
  s_train <- s[idx, ]; y_train <- y[idx]
  s_test <- s[-idx, ]; y_test <- y[-idx]

  pred_fun <- function(s_train, y_train, s_new) {
    fit <- lm(y_train ~ s_train[, 1] + s_train[, 2] +
                I(s_train[, 1]^2) + I(s_train[, 2]^2))
    cbind(1, s_new[, 1], s_new[, 2], s_new[, 1]^2, s_new[, 2]^2) %*% coef(fit)
  }

  out <- scp_geostatistical(s_train, y_train, s_test, pred_fun,
                             alpha = 0.1, seed = 7)

  expect_s3_class(out, "spconform")
  expect_equal(length(out$pred), nrow(s_test))
  expect_true(all(out$upper >= out$lower))

  cov <- coverage_report(out, y_test)
  # allow generous tolerance given small-sample / kernel approx
  expect_gt(cov$coverage, 0.70)
})

test_that("areal_neighbor_weights gives zero self-weight and decays with hops", {
  n <- 5
  adj <- matrix(0, n, n)
  for (i in 1:(n - 1)) { adj[i, i + 1] <- 1; adj[i + 1, i] <- 1 }
  w <- areal_neighbor_weights(1, adj)
  expect_equal(w[1], 0)
  expect_true(w[2] > w[3])
  expect_true(w[3] > w[4])
})

test_that("scp_areal returns valid intervals with default predictor", {
  set.seed(2)
  n <- 30
  adj <- matrix(0, n, n)
  for (i in 1:(n - 1)) { adj[i, i + 1] <- 1; adj[i + 1, i] <- 1 }
  y <- cumsum(rnorm(n)) + rnorm(n, sd = 0.2)

  out <- scp_areal(y, adjacency = adj, alpha = 0.2)
  expect_s3_class(out, "spconform")
  expect_equal(out$type, "areal")
  expect_true(all(out$upper >= out$lower))

  cov <- coverage_report(out, y)
  expect_gt(cov$coverage, 0.5)
})

test_that("input validation errors fire correctly", {
  expect_error(scp_geostatistical(matrix(1:4, 2), 1:3, matrix(1:2, 1), function(...) 1))
  expect_error(scp_areal(1:5, adjacency = matrix(0, 4, 4)))
})

test_that("S3 methods for spconform work correctly", {
  set.seed(123)
  s_tr <- matrix(runif(40), ncol = 2)
  y_tr <- rnorm(20)
  s_te <- matrix(runif(20), ncol = 2)
  y_te <- rnorm(10)

  pfun <- function(s_train, y_train, s_new) rep(mean(y_train), nrow(s_new))
  out <- scp_geostatistical(s_tr, y_tr, s_te, pfun, alpha = 0.1, seed = 1)

  # print & summary
  expect_output(print(out))
  expect_output(summary(out))

  # as.data.frame
  df <- as.data.frame(out)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), nrow(s_te))
  expect_true(all(c("pred", "lower", "upper", "width") %in% names(df)))

  # predict
  p_none <- predict(out, interval = "none")
  expect_equal(p_none, out$pred)
  p_mat <- predict(out, interval = "prediction")
  expect_true(is.matrix(p_mat))
  expect_equal(colnames(p_mat), c("fit", "lwr", "upr"))
  expect_equal(p_mat[, "fit"], out$pred)

  # residuals
  res_raw <- residuals(out, y_te, type = "response")
  expect_equal(res_raw, y_te - out$pred)
  res_abs <- residuals(out, y_te, type = "abs")
  expect_equal(res_abs, abs(y_te - out$pred))
  expect_error(residuals(out))

  # plot
  pdf(NULL)
  expect_silent(plot(out, y_true = y_te))
  dev.off()
})

test_that("diagnose suite returns classed object and supports plot/print", {
  set.seed(123)
  s_tr <- matrix(runif(40), ncol = 2)
  y_tr <- rnorm(20)
  s_te <- matrix(runif(20), ncol = 2)
  y_te <- rnorm(10)

  pfun <- function(s_train, y_train, s_new) rep(mean(y_train), nrow(s_new))
  out <- scp_geostatistical(s_tr, y_tr, s_te, pfun, alpha = 0.1, seed = 1)

  diag_res <- diagnose(out, y_te, s_te, plot = FALSE)
  expect_s3_class(diag_res, "spconform_diagnose")
  expect_true(!is.null(diag_res$marginal))
  expect_true(!is.null(diag_res$conditional))
  expect_true(!is.null(diag_res$boundary))
  expect_true(!is.null(diag_res$scores))

  # print method
  expect_output(print(diag_res))

  # plot method
  pdf(NULL)
  expect_silent(plot(diag_res))
  dev.off()
})



test_that("v0.1.2: single vector query point s0 works without error", {
  set.seed(42)
  s_tr <- matrix(runif(20), ncol = 2)
  y_tr <- rnorm(10)
  pfun <- function(s_train, y_train, s_new) rep(mean(y_train), nrow(s_new))

  # Single point query as a numeric vector (was a fatal dimension bug in v0.1.1)
  s0_vec <- c(0.45, 0.65)
  out_single <- scp_geostatistical(s_tr, y_tr, s0_vec, pfun, alpha = 0.1, seed = 1)
  expect_s3_class(out_single, "spconform")
  expect_equal(length(out_single$pred), 1)
  expect_true(out_single$upper >= out_single$lower)
})

test_that("v0.1.2: coverage_report computes Winkler Score and checks length mismatch", {
  set.seed(123)
  s_tr <- matrix(runif(40), ncol = 2)
  y_tr <- rnorm(20)
  s_te <- matrix(runif(20), ncol = 2)
  y_te <- rnorm(10)
  pfun <- function(s_train, y_train, s_new) rep(mean(y_train), nrow(s_new))
  out <- scp_geostatistical(s_tr, y_tr, s_te, pfun, alpha = 0.1, seed = 1)

  rep_res <- coverage_report(out, y_te)
  expect_true("winkler_score" %in% names(rep_res))
  expect_true(is.numeric(rep_res$winkler_score))
  expect_false(is.na(rep_res$winkler_score))

  # Length mismatch error
  expect_error(coverage_report(out, y_te[1:5]), "length")

  # Missing values handling
  y_te_na <- y_te
  y_te_na[1] <- NA
  rep_na_true <- coverage_report(out, y_te_na, na.rm = TRUE)
  expect_true(is.numeric(rep_na_true$coverage))
  expect_false(is.na(rep_na_true$coverage))

  rep_na_false <- coverage_report(out, y_te_na, na.rm = FALSE)
  expect_true(is.na(rep_na_false$coverage))
  expect_true(is.na(rep_na_false$winkler_score))
})

test_that("v0.1.2: k_neighbors parameter accelerates and functions correctly", {
  set.seed(99)
  n <- 50
  s_tr <- matrix(runif(2 * n), ncol = 2)
  y_tr <- rnorm(n)
  s_te <- matrix(runif(10), ncol = 2)
  pfun <- function(s_train, y_train, s_new) rep(mean(y_train), nrow(s_new))

  out_knn <- scp_geostatistical(s_tr, y_tr, s_te, pfun, alpha = 0.1, k_neighbors = 15, seed = 42)
  expect_s3_class(out_knn, "spconform")
  expect_equal(length(out_knn$pred), 5)
  expect_true(all(out_knn$upper >= out_knn$lower))
})

test_that("v0.1.2: defensive input validation detects NAs and bad arguments", {
  s_good <- matrix(1:4, 2, 2)
  y_good <- c(1, 2)
  pfun <- function(...) c(1, 1)

  # NAs in coordinates or response
  s_na <- s_good; s_na[1, 1] <- NA
  y_na <- y_good; y_na[1] <- NA
  expect_error(scp_geostatistical(s_na, y_good, s_good, pfun), "missing")
  expect_error(scp_geostatistical(s_good, y_na, s_good, pfun), "missing")
  expect_error(scp_geostatistical(s_good, y_good, s_na, pfun), "missing")

  # NAs in areal data
  adj <- matrix(c(0, 1, 1, 0), 2, 2)
  expect_error(scp_areal(c(1, NA), adjacency = adj), "missing")
  adj_na <- adj; adj_na[1, 2] <- NA
  expect_error(scp_areal(c(1, 2), adjacency = adj_na), "missing")
  expect_error(scp_areal(c(1, 2), X = matrix(c(1, NA), 2, 1), adjacency = adj), "missing")
})

test_that("v0.1.2: modular diagnostic functions work independently and programmatically", {
  set.seed(123)
  n <- 60
  s_tr <- matrix(runif(n * 2), n, 2)
  y_tr <- sin(3 * s_tr[, 1]) + rnorm(n, sd = 0.2)
  s_te <- matrix(runif(30 * 2), 30, 2)
  y_te <- sin(3 * s_te[, 1]) + rnorm(30, sd = 0.2)

  pfun <- function(s_train, y_train, s_new) rep(mean(y_train), nrow(s_new))
  fit <- scp_geostatistical(s_tr, y_tr, s_te, pfun, alpha = 0.1, seed = 1)

  # 1. spatial_moran
  m_res <- spatial_moran(fit, y_te, s_te)
  expect_s3_class(m_res, "spconform_moran")
  expect_true(is.numeric(m_res$moran_i))
  expect_true(is.numeric(m_res$p_value))
  expect_true(m_res$p_value >= 0 && m_res$p_value <= 1)

  # 2. boundary_diagnostics
  b_res <- boundary_diagnostics(fit, y_te, s_te)
  expect_s3_class(b_res, "spconform_boundary")
  expect_true(b_res$near_boundary$n > 0)
  expect_true(b_res$far_boundary$n > 0)
  expect_equal(b_res$near_boundary$n + b_res$far_boundary$n, 30)

  # 3. strata_diagnostics
  s_res <- strata_diagnostics(fit, y_te, s_te, n_bins = 2)
  expect_s3_class(s_res, "spconform_strata")
  expect_true(is.data.frame(s_res$strata))
  expect_true(nrow(s_res$strata) > 0)

  # 4. diagnose integration and as.data.frame
  d_res <- diagnose(fit, y_te, s_te, n_bins = 2, plot = FALSE)
  expect_s3_class(d_res, "spconform_diagnose")
  expect_s3_class(d_res$spatial_autocorr, "spconform_moran")
  expect_s3_class(d_res$boundary, "spconform_boundary")
  expect_s3_class(d_res$conditional, "spconform_strata")

  df_diag <- as.data.frame(d_res)
  expect_true(is.data.frame(df_diag))
  expect_true("coverage" %in% names(df_diag))
})
