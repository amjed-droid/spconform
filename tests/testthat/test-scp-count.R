test_that("scp_count constructs valid discrete integer intervals >= 0", {
  set.seed(42)
  n_tr <- 120
  s_tr <- matrix(runif(n_tr * 2), n_tr, 2)
  lambda_tr <- exp(1.2 + 1.5 * s_tr[, 1] - 0.8 * s_tr[, 2])
  y_tr <- stats::rpois(n_tr, lambda = lambda_tr)
  
  n_te <- 50
  s_te <- matrix(runif(n_te * 2), n_te, 2)
  lambda_te <- exp(1.2 + 1.5 * s_te[, 1] - 0.8 * s_te[, 2])
  y_te <- stats::rpois(n_te, lambda = lambda_te)
  
  pred_fun <- function(s_train, y_train, s_new) {
    df_tr <- data.frame(y = y_train, x = s_train[, 1], y_coord = s_train[, 2])
    df_nw <- data.frame(x = s_new[, 1], y_coord = s_new[, 2])
    fit <- stats::glm(y ~ x + y_coord, data = df_tr, family = stats::poisson)
    pmax(1e-4, as.numeric(stats::predict(fit, newdata = df_nw, type = "response")))
  }
  
  # 1. Pearson score (default)
  fit_cnt <- scp_count(s_tr, y_tr, s_te, pred_fun, alpha = 0.10, score_type = "pearson")
  expect_s3_class(fit_cnt, "spconform_count")
  expect_s3_class(fit_cnt, "spconform")
  
  # Strict non-negativity and integer checks
  expect_true(all(fit_cnt$lower >= 0L))
  expect_true(all(fit_cnt$upper >= fit_cnt$lower))
  expect_true(is.integer(fit_cnt$lower))
  expect_true(is.integer(fit_cnt$upper))
  
  # Coverage evaluation
  rep_cnt <- coverage_report(fit_cnt, y_te)
  expect_true(rep_cnt$coverage >= 0.85) # Targets 90%
  
  # 2. Anscombe score
  fit_ans <- scp_count(s_tr, y_tr, s_te, pred_fun, alpha = 0.10, score_type = "anscombe")
  expect_true(all(fit_ans$lower >= 0L))
  rep_ans <- coverage_report(fit_ans, y_te)
  expect_true(rep_ans$coverage >= 0.85)
})

test_that("scp_count handles overdispersion, S3 methods, and input errors", {
  set.seed(123)
  n_tr <- 100
  s_tr <- matrix(runif(n_tr * 2), n_tr, 2)
  # Overdispersed negative binomial-like counts
  y_tr <- stats::rnbinom(n_tr, size = 2, mu = exp(1.5 + s_tr[, 1]))
  
  s_te <- matrix(runif(30 * 2), 30, 2)
  
  pred_fun <- function(s_train, y_train, s_new) {
    df_tr <- data.frame(y = y_train, x = s_train[, 1], y_coord = s_train[, 2])
    df_nw <- data.frame(x = s_new[, 1], y_coord = s_new[, 2])
    fit <- stats::glm(y ~ x + y_coord, data = df_tr, family = stats::quasipoisson)
    pmax(1e-4, as.numeric(stats::predict(fit, newdata = df_nw, type = "response")))
  }
  
  fit_nb <- scp_count(s_tr, y_tr, s_te, pred_fun, alpha = 0.10)
  expect_true(fit_nb$dispersion >= 0)
  
  # Test coverage_report for count modality
  y_te <- stats::rnbinom(30, size = 2, mu = exp(1.5 + s_te[, 1]))
  rep_cnt <- coverage_report(fit_nb, y_te)
  expect_s3_class(rep_cnt, "spconform_coverage_report")
  expect_true(rep_cnt$coverage >= 0.70)
  expect_true(rep_cnt$winkler_score > 0)
  expect_true(rep_cnt$mean_span_width >= 1)
  expect_true(rep_cnt$zero_lower_rate >= 0 && rep_cnt$zero_lower_rate <= 1)
  expect_output(print(rep_cnt), "spconform Coverage & Efficiency Audit")

  # Test S3 print and summary
  expect_output(print(fit_nb))
  expect_output(summary(fit_nb))
  
  # Test S3 plot
  temp_pdf <- tempfile(fileext = ".pdf")
  pdf(temp_pdf)
  expect_silent(plot(fit_nb))
  dev.off()
  unlink(temp_pdf)
  
  # Input errors
  expect_error(scp_count(s_tr, c(-1, y_tr[-1]), s_te, pred_fun)) # Negative counts
  expect_error(scp_count(s_tr, y_tr, s_te, pred_fun, alpha = -0.05)) # Bad alpha
  expect_error(scp_count(s_tr, y_tr[1:10], s_te, pred_fun)) # Length mismatch
})
