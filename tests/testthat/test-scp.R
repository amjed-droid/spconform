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

