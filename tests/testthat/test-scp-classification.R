test_that("scp_classification works correctly on binary spatial data", {
  set.seed(42)
  n_tr <- 150
  s_tr <- matrix(runif(n_tr * 2), n_tr, 2)
  prob_healthy <- 1 / (1 + exp(-20 * (s_tr[, 1] + s_tr[, 2] - 1.0)))
  y_tr <- factor(ifelse(runif(n_tr) < prob_healthy, "Healthy", "Diseased"))
  
  n_te <- 50
  s_te <- matrix(runif(n_te * 2), n_te, 2)
  prob_te <- 1 / (1 + exp(-20 * (s_te[, 1] + s_te[, 2] - 1.0)))
  y_te <- factor(ifelse(runif(n_te) < prob_te, "Healthy", "Diseased"))
  
  pred_fun <- function(s_train, y_train, s_new) {
    df_tr <- data.frame(y = y_train, x = s_train[, 1], y_coord = s_train[, 2])
    df_nw <- data.frame(x = s_new[, 1], y_coord = s_new[, 2])
    fit <- stats::glm(y ~ x + y_coord, data = df_tr, family = stats::binomial)
    p_healthy <- stats::predict(fit, newdata = df_nw, type = "response")
    cbind(Healthy = p_healthy, Diseased = 1 - p_healthy)
  }
  
  # 1. Fit LAC classification (Least Likely / Top-class)
  fit_lac <- scp_classification(s_tr, y_tr, s_te, pred_fun, alpha = 0.10, method = "LAC")
  expect_s3_class(fit_lac, "spconform_classification")
  rep_lac <- coverage_report(fit_lac, y_te)
  expect_true(rep_lac$coverage >= 0.85)
  expect_true(rep_lac$singleton_rate >= 0.50)
  
  # 2. Fit APS classification (Adaptive Prediction Sets)
  fit_aps <- scp_classification(s_tr, y_tr, s_te, pred_fun, alpha = 0.10, method = "APS")
  expect_s3_class(fit_aps, "spconform_classification")
  rep_aps <- coverage_report(fit_aps, y_te)
  expect_true(rep_aps$coverage >= 0.85)
  expect_true(rep_aps$avg_set_size >= 1.0)
})

test_that("scp_classification works on multiclass spatial data and S3 methods work", {
  set.seed(123)
  n_tr <- 150
  s_tr <- matrix(runif(n_tr * 2), n_tr, 2)
  zone <- ifelse(s_tr[, 1] < 0.33, "Zone_A",
          ifelse(s_tr[, 1] < 0.66, "Zone_B", "Zone_C"))
  y_tr <- factor(zone)
  
  n_te <- 40
  s_te <- matrix(runif(n_te * 2), n_te, 2)
  zone_te <- ifelse(s_te[, 1] < 0.33, "Zone_A",
             ifelse(s_te[, 1] < 0.66, "Zone_B", "Zone_C"))
  y_te <- factor(zone_te)
  
  pred_fun_multi <- function(s_train, y_train, s_new) {
    cls <- levels(y_train)
    p_mat <- matrix(0.1, nrow = nrow(s_new), ncol = length(cls))
    colnames(p_mat) <- cls
    for (i in 1:nrow(s_new)) {
      if (s_new[i, 1] < 0.33) p_mat[i, "Zone_A"] <- 0.8
      else if (s_new[i, 1] < 0.66) p_mat[i, "Zone_B"] <- 0.8
      else p_mat[i, "Zone_C"] <- 0.8
    }
    p_mat / rowSums(p_mat)
  }
  
  fit_multi <- scp_classification(s_tr, y_tr, s_te, pred_fun_multi, alpha = 0.10)
  expect_equal(fit_multi$classes, c("Zone_A", "Zone_B", "Zone_C"))
  
  # Test S3 print and summary
  expect_output(print(fit_multi))
  expect_output(summary(fit_multi))
  
  # Test S3 plot
  temp_pdf <- tempfile(fileext = ".pdf")
  pdf(temp_pdf)
  expect_silent(plot(fit_multi))
  dev.off()
  unlink(temp_pdf)
  
  # Input errors
  expect_error(scp_classification(s_tr, y_tr, s_te, pred_fun_multi, alpha = 1.5))
  expect_error(scp_classification(s_tr, y_tr[1:10], s_te, pred_fun_multi))
})
