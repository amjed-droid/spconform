test_that("scp_multivariate input validation and errors", {
  set.seed(42)
  n <- 30
  s <- matrix(runif(n * 2), n, 2)
  y_1d <- matrix(rnorm(n), n, 1)

  pfun <- function(st, yt, sn) cbind(rep(0, nrow(sn)), rep(0, nrow(sn)))

  # d < 2 error
  expect_error(scp_multivariate(s, y_1d, s[1:5, ], pfun), "requires at least 2 response dimensions")

  # Mismatched rows error
  expect_error(scp_multivariate(s, matrix(rnorm(20 * 2), 20, 2), s[1:5, ], pfun), "nrow\\(y_train\\) must equal nrow\\(s_train\\)")

  # Invalid alpha
  expect_error(scp_multivariate(s, matrix(rnorm(n * 2), n, 2), s[1:5, ], pfun, alpha = 1.2), "alpha must be strictly between 0 and 1")
})

test_that("scp_multivariate generates valid bivariate ellipsoids and bounding intervals", {
  set.seed(123)
  n <- 80
  s <- matrix(runif(n * 2, 0, 10), n, 2)
  y1 <- sin(s[, 1]) + rnorm(n, sd = 0.2)
  y2 <- cos(s[, 2]) + 0.6 * y1 + rnorm(n, sd = 0.2)
  y <- cbind(zinc = y1, lead = y2)

  s_tr <- s[1:60, ]
  y_tr <- y[1:60, ]
  s_te <- s[61:80, ]
  y_te <- y[61:80, ]

  pfun <- function(st, yt, sn) {
    df_fit <- data.frame(y1 = yt[, 1], y2 = yt[, 2], x = st[, 1], y = st[, 2])
    df_nw  <- data.frame(x = sn[, 1], y = sn[, 2])
    m1 <- stats::lm(y1 ~ x + y, data = df_fit)
    m2 <- stats::lm(y2 ~ x + y, data = df_fit)
    cbind(as.numeric(stats::predict(m1, newdata = df_nw)),
          as.numeric(stats::predict(m2, newdata = df_nw)))
  }

  fit <- scp_multivariate(s_tr, y_tr, s_te, pfun, alpha = 0.10)

  expect_s3_class(fit, "spconform_multivariate")
  expect_s3_class(fit, "spconform")
  expect_equal(fit$dimension, 2)
  expect_equal(nrow(fit$pred), 20)
  expect_equal(ncol(fit$pred), 2)
  expect_equal(nrow(fit$lower), 20)
  expect_equal(nrow(fit$upper), 20)

  # Projected bounding intervals satisfy lower <= upper
  expect_true(all(fit$upper >= fit$lower))

  # Ellipsoid volumes and radii are positive
  expect_true(all(fit$quantiles > 0))
  expect_true(all(fit$ellipsoid_volumes > 0))

  # Covariance matrices are 2x2 positive definite
  expect_equal(length(fit$cov_matrices), 20)
  for (cov_m in fit$cov_matrices) {
    expect_equal(dim(cov_m), c(2, 2))
    expect_true(det(cov_m) > 0)
  }

  # Coverage report evaluation
  cov_rep <- multivariate_coverage_report(fit, y_te)
  expect_s3_class(cov_rep, "spconform_multicov_report")
  expect_equal(cov_rep$nominal_confidence, 0.90)
  expect_true(cov_rep$joint_coverage >= 0.70) # Empirical finite-sample coverage should be high
  expect_true(all(cov_rep$marginal_coverage >= cov_rep$joint_coverage)) # Marginals must be >= joint

  # Print methods execute cleanly
  expect_output(print(fit), "=== Spatially Localized Multivariate Conformal Prediction ===")
  expect_output(print(cov_rep), "=== Multivariate Conformal Coverage Evaluation ===")
})

test_that("scp_multivariate scales to d = 3 (trivariate) with different kernels", {
  set.seed(456)
  n <- 90
  s <- matrix(runif(n * 2, 0, 5), n, 2)
  y1 <- s[, 1] + rnorm(n, sd = 0.3)
  y2 <- s[, 2] + 0.3 * y1 + rnorm(n, sd = 0.3)
  y3 <- 0.5 * s[, 1] - 0.2 * s[, 2] + rnorm(n, sd = 0.3)
  y <- cbind(cadmium = y1, copper = y2, lead = y3)

  pfun3 <- function(st, yt, sn) {
    df_fit <- data.frame(y1 = yt[, 1], y2 = yt[, 2], y3 = yt[, 3], x = st[, 1], y = st[, 2])
    df_new <- data.frame(x = sn[, 1], y = sn[, 2])
    m1 <- stats::lm(y1 ~ x + y, data = df_fit)
    m2 <- stats::lm(y2 ~ x + y, data = df_fit)
    m3 <- stats::lm(y3 ~ x + y, data = df_fit)
    cbind(stats::predict(m1, newdata = df_new),
          stats::predict(m2, newdata = df_new),
          stats::predict(m3, newdata = df_new))
  }

  fit_exp <- scp_multivariate(s[1:65, ], y[1:65, ], s[66:90, ], pfun3,
                              alpha = 0.15, kernel = "exponential", h = 2.0)
  expect_equal(fit_exp$dimension, 3)
  expect_true(all(fit_exp$upper >= fit_exp$lower))

  fit_bi <- scp_multivariate(s[1:65, ], y[1:65, ], s[66:90, ], pfun3,
                             alpha = 0.15, kernel = "bisquare", h = 3.5)
  expect_equal(fit_bi$dimension, 3)
  expect_true(all(fit_bi$upper >= fit_bi$lower))
})
