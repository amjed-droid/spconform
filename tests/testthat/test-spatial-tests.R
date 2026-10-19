test_that("spatial_conditional_coverage works as expected under H0 and H1", {
  set.seed(42)
  n <- 80
  coords <- matrix(runif(n * 2), n, 2)
  
  mock_fit <- list(
    pred = rep(0, n),
    lower = rep(-1, n),
    upper = rep(1, n),
    alpha = 0.1,
    type = "geostatistical",
    coords = coords
  )
  class(mock_fit) <- "spconform"
  
  # 1. Uniform random coverage (H0 holds)
  y_null <- rnorm(n, mean = 0, sd = 0.6)
  test_null <- spatial_conditional_coverage(mock_fit, y_null, coords, test = "LRT")
  expect_s3_class(test_null, "spconform_coverage_test")
  expect_true(test_null$p_value > 0.05)
  expect_match(test_null$decision, "Spatial conditional coverage holds")
  
  # 2. Geographically biased coverage (H1 holds)
  prob_miss <- ifelse(coords[, 1] > 0.6, 0.8, 0.1)
  is_miss <- rbinom(n, size = 1, prob = prob_miss)
  y_biased <- ifelse(is_miss == 1, 5.0, 0.0)
  test_biased <- spatial_conditional_coverage(mock_fit, y_biased, coords, degree = 1, test = "LRT")
  expect_true(test_biased$p_value < 0.05)
  expect_match(test_biased$decision, "Spatial conditional coverage violated")
  
  # 3. Wald test
  test_wald <- spatial_conditional_coverage(mock_fit, y_null, coords, test = "Wald")
  expect_true(test_wald$df %in% c(2, 5))
  
  # 4. Homogeneous perfect coverage edge case
  y_perfect <- rep(0, n)
  test_perfect <- spatial_conditional_coverage(mock_fit, y_perfect, coords)
  expect_equal(test_perfect$p_value, 1.0)
  expect_equal(test_perfect$statistic, 0)
})

test_that("spatial_conditional_coverage handles input errors correctly", {
  mock_fit <- list(pred = rep(0, 10), lower = rep(-1, 10), upper = rep(1, 10), alpha = 0.1)
  class(mock_fit) <- "spconform"
  
  expect_error(spatial_conditional_coverage(mock_fit, rep(0, 5), matrix(0, 10, 2)))
  expect_error(spatial_conditional_coverage(mock_fit, rep(0, 10), NULL))
  expect_error(spatial_conditional_coverage(list(), rep(0, 10), matrix(0, 10, 2)))
})

test_that("spatial_geary computes correct statistics under H0 and H1", {
  set.seed(123)
  n <- 100
  coords <- matrix(runif(n * 2), n, 2)
  
  # 1. Independent random noise (H0: C ~ 1)
  x_iid <- rnorm(n)
  geary_iid <- spatial_geary(x_iid, coords)
  expect_s3_class(geary_iid, "spconform_geary")
  expect_equal(geary_iid$expected, 1.0)
  expect_true(geary_iid$variance > 0)
  expect_true(abs(geary_iid$C - 1.0) < 0.15)
  expect_true(geary_iid$p_value > 0.05)
  
  # 2. Strong spatial clustering / smooth gradient (H1: C < 1, positive spatial autocorrelation)
  x_clustered <- coords[, 1] * 4 + coords[, 2] * 4 + rnorm(n, sd = 0.1)
  geary_clustered <- spatial_geary(x_clustered, coords, alternative = "greater")
  expect_true(geary_clustered$C < 0.75)
  expect_true(geary_clustered$z_score > 5.0)
  expect_true(geary_clustered$p_value < 0.0001)
  
  # 3. Input validation
  expect_error(spatial_geary(c(1, 2), coords[1:2, ]))
  expect_error(spatial_geary(rep(5, 10), coords[1:10, ]))
})

test_that("spatial_lisa computes correct statistics and classifies quadrants", {
  set.seed(42)
  n <- 80
  coords <- matrix(runif(n * 2), n, 2)
  
  # Hotspot in the Northeast quadrant (x > 0.5 & y > 0.5)
  x_val <- ifelse(coords[, 1] > 0.5 & coords[, 2] > 0.5, 4.0, 0.0) + rnorm(n, sd = 0.3)
  
  lisa_res <- spatial_lisa(x_val, coords, p_adjust_method = "fdr")
  expect_s3_class(lisa_res, "spconform_lisa")
  expect_equal(length(lisa_res$Ii), n)
  expect_equal(length(lisa_res$quadrant), n)
  expect_true(all(lisa_res$variance > 0))
  
  # Hotspot detection: High-High clusters must be detected in NE
  counts <- table(lisa_res$quadrant)
  expect_true(counts["High-High"] > 0)
  
  # Check Bonferroni vs FDR
  lisa_bonf <- spatial_lisa(x_val, coords, p_adjust_method = "bonferroni")
  expect_true(all(lisa_bonf$p_adjusted >= lisa_res$p_adjusted))
  
  # Test plot method does not error
  temp_pdf <- tempfile(fileext = ".pdf")
  pdf(temp_pdf)
  expect_silent(plot(lisa_res))
  dev.off()
  unlink(temp_pdf)
  
  # Input errors
  expect_error(spatial_lisa(1:3, coords[1:3, ])) # n < 5
  expect_error(spatial_lisa(rep(1, 10), coords[1:10, ])) # zero variance
})


test_that("spatial_anisotropy works correctly under isotropy and anisotropy", {
  set.seed(42)
  n <- 80
  coords <- matrix(runif(n * 2, 0, 10), n, 2)
  
  # 1. Isotropic noise (H0 holds)
  x_iso <- rnorm(n)
  aniso_iso <- spatial_anisotropy(x_iso, coords = coords, n_perm = 99)
  expect_s3_class(aniso_iso, "spconform_anisotropy")
  expect_true(aniso_iso$p_value > 0.05)
  expect_match(aniso_iso$decision, "Spatial isotropy confirmed")
  
  # 2. Strong directional anisotropy (strong variation along X axis, constant along Y)
  x_aniso <- coords[, 1] * 3 + rnorm(n, sd = 0.2)
  aniso_dir <- spatial_anisotropy(x_aniso, coords = coords, n_perm = 99)
  expect_true(aniso_dir$anisotropy_ratio > 1.5)
  expect_true(aniso_dir$p_value < 0.05)
  expect_match(aniso_dir$decision, "Directional anisotropy detected")
  
  # 3. Test plot method executes silently
  temp_pdf <- tempfile(fileext = ".pdf")
  pdf(temp_pdf)
  expect_silent(plot(aniso_dir))
  dev.off()
  unlink(temp_pdf)
  
  # 4. Input error
  expect_error(spatial_anisotropy(1:4, coords = coords[1:4, ])) # n < 8
})
