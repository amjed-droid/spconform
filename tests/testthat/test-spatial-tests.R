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


test_that("spatial_mantel correctly evaluates matrix correlation and whitening", {
  set.seed(42)
  n <- 50
  coords <- matrix(runif(n * 2, 0, 10), n, 2)

  # 1. Spatially independent noise (H0 holds: spatial whitening)
  scores_indep <- rnorm(n)
  res_indep <- spatial_mantel(scores_indep, coords = coords, n_permutations = 199, seed = 123)
  expect_s3_class(res_indep, "spconform_mantel")
  expect_equal(res_indep$n, n)
  expect_equal(res_indep$n_pairs, n * (n - 1) / 2)
  expect_true(res_indep$p_value >= 0.05)
  expect_true(res_indep$is_whitened)

  # 2. Spatially correlated scores (H0 rejected: clustering)
  scores_corr <- sin(coords[, 1] * 0.4) + cos(coords[, 2] * 0.4) + rnorm(n, sd = 0.1)
  res_corr <- spatial_mantel(scores_corr, coords = coords, n_permutations = 199, seed = 123)
  expect_true(res_corr$mantel_r > 0.20)
  expect_true(res_corr$p_value < 0.05)
  expect_false(res_corr$is_whitened)

  # 3. Integration with spconform object
  pfun <- function(st, yt, sn) rep(mean(yt), nrow(sn))
  fit_sp <- scp_geostatistical(coords[1:35, ], scores_indep[1:35], coords[36:50, ], pfun, alpha = 0.1)
  res_obj <- spatial_mantel(fit_sp, y_true = scores_indep[36:50], n_permutations = 99, seed = 42)
  expect_s3_class(res_obj, "spconform_mantel")

  # 4. Print and Plot methods execute cleanly
  expect_output(print(res_indep), "spconform: Spatial Mantel Test")

  temp_pdf <- tempfile(fileext = ".pdf")
  pdf(temp_pdf)
  expect_silent(plot(res_indep))
  dev.off()
  unlink(temp_pdf)

  # 5. Error handling
  expect_error(spatial_mantel(rnorm(3), coords = coords[1:3, ]), "Need at least 5 spatial locations")
  expect_error(spatial_mantel(rnorm(10), coords = coords[1:5, ]), "nrow\\(coords\\) must equal length")
})


test_that("spatial_cross_l evaluates bivariate hotspot attraction and repulsion", {
  set.seed(42)
  n <- 50
  coords <- matrix(runif(n * 2, 0, 10), n, 2)
  y_true <- sin(coords[, 1]) + rnorm(n, sd = 0.3)
  pfun <- function(st, yt, sn) rep(mean(yt), nrow(sn))

  # 1. Standard execution with spconform object
  fit <- scp_geostatistical(coords[1:35, ], y_true[1:35], coords[36:50, ], pfun, alpha = 0.1)
  res_l <- spatial_cross_l(fit, y_true = y_true[36:50], n_sim = 49, seed = 123)

  expect_s3_class(res_l, "spconform_cross_l")
  expect_equal(res_l$n_total, 15)
  if (res_l$n_miss > 0) {
    expect_true(length(res_l$r) > 0)
    expect_true(length(res_l$l_obs) == length(res_l$r))
    expect_true(all(res_l$l_upper >= res_l$l_lower))
    expect_true(res_l$pattern %in% c("Spatial Attraction (Honest Uncertainty)",
                                     "Spatial Repulsion (Overconfidence)",
                                     "Spatial Randomness"))
  } else {
    expect_equal(res_l$pattern, "Perfect Empirical Coverage")
  }

  # 2. Test perfect coverage edge case
  bounds_perfect <- cbind(y_true[36:50] - 10, y_true[36:50] + 10)
  res_perf <- spatial_cross_l(bounds_perfect, coords = coords[36:50, ], y_true = y_true[36:50])
  expect_equal(res_perf$n_miss, 0)
  expect_equal(res_perf$pattern, "Perfect Empirical Coverage")

  # 3. Print and plot methods execute cleanly
  expect_output(print(res_l), "spconform: Bivariate Cross-L Ripley Hotspot Audit")

  temp_pdf <- tempfile(fileext = ".pdf")
  pdf(temp_pdf)
  expect_silent(plot(res_l))
  expect_silent(plot(res_perf))
  dev.off()
  unlink(temp_pdf)

  # 4. Input validation
  expect_error(spatial_cross_l(fit), "'y_true' must be provided")
  expect_error(spatial_cross_l(bounds_perfect[1:5, ], coords = coords[1:5, ], y_true = y_true[1:5]),
               "Need at least 10 spatial locations")
})

test_that("conformal_ce_surface computes local efficiency-coverage trade-off correctly", {
  set.seed(42)
  n <- 50
  coords <- matrix(runif(n * 2, 0, 10), n, 2)
  y_true <- sin(coords[, 1]) + rnorm(n, sd = 0.3)
  pfun <- function(st, yt, sn) rep(mean(yt), nrow(sn))

  # 1. Geostatistical continuous evaluation
  fit <- scp_geostatistical(coords[1:35, ], y_true[1:35], coords[36:50, ], pfun, alpha = 0.1)
  ce <- conformal_ce_surface(fit, y_true = y_true[36:50], grid_resolution = 20)

  expect_s3_class(ce, "spconform_ce_surface")
  expect_equal(ce$n, 15)
  expect_equal(length(ce$local_coverage), 15)
  expect_equal(length(ce$local_width), 15)
  expect_equal(levels(ce$zones), c("Optimal", "Conservative", "Vulnerable", "Deficit"))
  expect_equal(sum(ce$zone_distribution$Percentage), 100, tolerance = 0.5)
  expect_true(is.numeric(ce$tradeoff_cor))
  expect_true(is.numeric(ce$ecdi))
  expect_true(!is.null(ce$grid_surface))
  expect_equal(dim(ce$grid_surface$coverage_matrix), c(20, 20))

  # Alias check
  ce_alias <- spatial_ce_surface(fit, y_true = y_true[36:50], grid_resolution = 0)
  expect_s3_class(ce_alias, "spconform_ce_surface")
  expect_null(ce_alias$grid_surface)

  # 2. Print and plot methods execute cleanly across all types
  expect_output(print(ce), "spconform: Local Efficiency-Coverage Trade-off Surface")

  temp_pdf <- tempfile(fileext = ".pdf")
  pdf(temp_pdf)
  expect_silent(plot(ce, type = "phase"))
  expect_silent(plot(ce, type = "map"))
  expect_silent(plot(ce, type = "coverage"))
  expect_silent(plot(ce, type = "width"))
  expect_silent(plot(ce_alias, type = "coverage"))
  expect_silent(plot(ce_alias, type = "width"))
  dev.off()
  unlink(temp_pdf)

  # 3. Classification modality
  fit_cl <- list(
    prediction_sets = list(c("A"), c("A", "B"), c("B"), c("A"), c("B")),
    set_sizes = c(1, 2, 1, 1, 1),
    alpha = 0.1,
    s0 = coords[1:5, ]
  )
  class(fit_cl) <- "spconform_classification"
  ce_cl <- conformal_ce_surface(fit_cl, y_true = c("A", "A", "B", "B", "A"), grid_resolution = 10)
  expect_s3_class(ce_cl, "spconform_ce_surface")
  expect_equal(ce_cl$efficiency_name, "Prediction Set Size")

  # 4. Count modality
  fit_cnt <- list(
    pred = rep(5, 5),
    lower = c(2, 3, 1, 4, 2),
    upper = c(8, 7, 6, 9, 7),
    alpha = 0.1,
    s0 = coords[1:5, ]
  )
  class(fit_cnt) <- "spconform_count"
  ce_cnt <- conformal_ce_surface(fit_cnt, y_true = c(4, 5, 2, 12, 6), grid_resolution = 10)
  expect_s3_class(ce_cnt, "spconform_ce_surface")
  expect_equal(ce_cnt$efficiency_name, "Count Span Width")

  # 5. Raw matrix of bounds
  mat_bounds <- cbind(lower = c(-1, -2, -1, -3, -2), upper = c(1, 2, 1, 3, 2))
  ce_mat <- conformal_ce_surface(mat_bounds, y_true = c(0, 0, 0, 5, 0), coords = coords[1:5, ])
  expect_s3_class(ce_mat, "spconform_ce_surface")

  # 6. Error handling
  expect_error(conformal_ce_surface(mat_bounds[1:3, ], y_true = c(0, 0, 0), coords = coords[1:3, ]),
               "Need at least 5 spatial locations")
  expect_error(conformal_ce_surface(mat_bounds, y_true = c(0, 0), coords = coords[1:5, ]),
               "Dimensions of coords, y_true, and prediction intervals/sets must match")
  expect_error(conformal_ce_surface(list(), y_true = c(0, 0)),
               "must be an 'spconform' object")
})

