statistical_blocks = function() {
  t = seq_len(24)
  x = cbind(t / 8, sin(t), cos(t / 3))
  y = cbind(2 * x[, 1] + cos(t), x[, 2] - x[, 3], sin(t / 2))
  list(x = x, y = y)
}

test_that("PLS reconstructs between-block covariance and centered scores", {
  blocks = statistical_blocks()
  x = blocks$x
  y = blocks$y
  fit = pls(as.data.frame(x), y, perm = 0)
  expect_s3_class(fit, "pls_fit")
  expect_equal(fit$U %*% diag(fit$D) %*% t(fit$V), cov(x, y))
  expect_equal(crossprod(fit$U), diag(3))
  expect_equal(crossprod(fit$V), diag(3))
  expect_equal(unname(fit$XScores %*% t(fit$U)), unname(sweep(x, 2, colMeans(x))))
  expect_equal(unname(fit$YScores %*% t(fit$V)), unname(sweep(y, 2, colMeans(y))))
  expect_equal(sum(fit$D^2), sum(cov(x, y)^2))
  expect_equal(sum(fit$percentage_squared_covariance), 100)
  expect_equal(fit$CorrXScoresYScores, diag(cor(fit$XScores, fit$YScores)))
  translated = pls(sweep(x, 2, c(5, -2, 8), "+"), y + 10, perm = 0)
  expect_equal(translated$D, fit$D)
  expect_equal(unname(fit$OriginalData$x), unname(x))
  expect_equal(unname(fit$x_center), unname(colMeans(x)))
  expect_output(print(fit), "Two-block PLS analysis")
  expect_output(summary(fit), "percentage_squared_covariance")
})

test_that("PLS supports unequal block widths and rejects incompatible inputs", {
  blocks = statistical_blocks()
  fit = pls(blocks$x[, 1:2], blocks$y, perm = 0)
  expect_equal(dim(fit$U), c(2L, 2L))
  expect_equal(dim(fit$V), c(3L, 2L))
  expect_equal(fit$U %*% diag(fit$D) %*% t(fit$V), cov(blocks$x[, 1:2], blocks$y))
  expect_error(pls(blocks$x[-1, ], blocks$y, perm = 0), "same number of rows")
  expect_error(pls(blocks$x[, FALSE, drop = FALSE], blocks$y, perm = 0), "at least one column")
})

test_that("PLS permutations preserve fits and give reproducible attainable p values", {
  skip_if_not_installed("future.apply")
  skip_if_not_installed("future")
  withr::local_seed(924)
  old_plan = future::plan(future::sequential)
  withr::defer(future::plan(old_plan))
  blocks = statistical_blocks()
  fit = pls(blocks$x, blocks$y, perm = 40)
  set.seed(924)
  repeated = pls(blocks$x, blocks$y, perm = 40)
  expect_equal(fit, repeated)
  expect_equal(fit$D, pls(blocks$x, blocks$y, perm = 0)$D)
  expect_equal(unname(fit$global_significance_RV["observed_RV"]), EscoufierRV(blocks$x, blocks$y))
  p = c(fit$global_significance_RV["RV_p_value"],
        fit$singular_axis_significance$pvalue_singular_values,
        fit$singular_axis_significance$pvalue_correlation_PLS_scores)
  expect_true(all(p >= 1 / 41 & p <= 1))
  expect_equal(unname(p * 41), unname(round(p * 41)))
  expect_output(print(fit), "Observed RV=")
  expect_output(summary(fit), "Significance of each singular axis")
  axis_only = pls(blocks$x, blocks$y, perm = 40, global_RV_test = FALSE)
  expect_false("global_significance_RV" %in% names(axis_only))
  expect_output(print(axis_only), "correlation_between_scores")
  expect_output(summary(axis_only), "Significance of each singular axis")
})

test_that("major-axis predictions reconstruct exactly related blocks without sign assumptions", {
  x = statistical_blocks()$x
  y = sweep(2 * x, 2, c(3, -2, 1), "+")
  fit = pls(x, y, perm = 0)
  new_x = x[c(2, 9, 17), , drop = FALSE] + 0.3
  new_y = sweep(2 * new_x, 2, c(3, -2, 1), "+")
  for (scaled in c(TRUE, FALSE)) {
    prediction = pls_major_axis(fit, new_x, new_y, axes_to_use = 3, scale_PLS = scaled)
    original = prediction$original_major_axis_predictions_reversed
    expect_equal(unname(original$Block1), unname(x))
    expect_equal(unname(original$Block2), unname(y))
    projected = prediction$new_data_results
    expect_equal(unname(projected$new_data_Block1_proj_prediction_revert), unname(new_x))
    expect_equal(unname(projected$new_data_Block2_proj_prediction_revert), unname(new_y))
    expect_length(prediction$original_major_axis_projection, 3)
  }
  first_axis = pls_major_axis(fit, axes_to_use = 1)
  expected = sweep(fit$XScores[, 1, drop = FALSE] %*% t(fit$U[, 1, drop = FALSE]),
                   2, colMeans(x), "+")
  expect_equal(unname(first_axis$original_major_axis_predictions_reversed$Block1), unname(expected))
  expect_error(pls_major_axis(list()), "class 'pls_fit'")
  expect_error(pls_major_axis(fit, axes_to_use = 0), "axes_to_use")
  expect_error(pls_major_axis(fit, axes_to_use = 4), "axes_to_use")
  expect_error(pls_major_axis(fit, new_x[, 1:2], new_y), "number of variables")
})

test_that("major-axis predictions remove off-axis variation by geometric projection", {
  x = matrix(seq(-4, 4), ncol = 1)
  y = matrix(2 * x[, 1] + c(2, -1, 0, 1, -2, 1, 0, -1, 2) + 5, ncol = 1)
  fit = pls(x, y, perm = 0)
  data = cbind(x, y)
  center = colMeans(data)
  centered = sweep(data, 2, center)
  new_data = rbind(c(-2, 5), c(1, 4), c(3, 7))
  new_centered = sweep(new_data, 2, center)
  covariance = cov(data)
  for (scaled in c(TRUE, FALSE)) {
    if (scaled) {
      standard_deviations = sqrt(diag(covariance))
      standardized = sweep(centered, 2, standard_deviations, "/")
      new_standardized = sweep(new_centered, 2, standard_deviations, "/")
      # Positive correlation and equal standardized variances give a 45-degree major axis.
      expected_centered = outer(rowMeans(standardized), standard_deviations)
      expected_new_centered = outer(rowMeans(new_standardized), standard_deviations)
    } else {
      angle = atan2(2 * covariance[1, 2], covariance[1, 1] - covariance[2, 2]) / 2
      direction = c(cos(angle), sin(angle)) # Closed-form principal direction of a 2D ellipse.
      expected_centered = centered %*% tcrossprod(direction)
      expected_new_centered = new_centered %*% tcrossprod(direction)
    }
    expected = sweep(expected_centered, 2, center, "+")
    expected_new = sweep(expected_new_centered, 2, center, "+")
    prediction = pls_major_axis(fit, new_data[, 1, drop = FALSE],
                                new_data[, 2, drop = FALSE], scale_PLS = scaled)
    actual = prediction$original_major_axis_predictions_reversed
    actual_new = prediction$new_data_results
    expect_equal(unname(cbind(actual$Block1, actual$Block2)), unname(expected))
    expect_equal(unname(cbind(actual_new$new_data_Block1_proj_prediction_revert,
                              actual_new$new_data_Block2_proj_prediction_revert)),
                 unname(expected_new))
    expect_gt(sum((expected - data)^2), 0.1)
    expect_gt(sum((expected_new - new_data)^2), 0.1)
  }
})

test_that("paired Hotelling statistics agree with a multivariate intercept Wilks test", {
  blocks = statistical_blocks()
  a = blocks$x[, 1:2]
  b = a + blocks$y[, 1:2] / 3
  rownames(a) = rownames(b) = paste0("specimen_", seq_len(nrow(a)))
  expect_warning({ result = repeated_measures_test(a, b) }, "used for matching")
  difference = b - a
  residuals = residuals(lm(difference ~ 1))
  lambda = det(crossprod(residuals)) / det(crossprod(difference)) # Wilks' lambda for a mean of zero.
  n = nrow(a)
  p = ncol(a)
  expected_f = (n - p) / p * (1 - lambda) / lambda
  expect_equal(unname(result["Fstat"]), expected_f)
  expect_equal(unname(result["HotellingT2"]), expected_f * p * (n - 1) / (n - p))
  expect_equal(unname(result["p_value"]), pf(expected_f, p, n - p, lower.tail = FALSE))
  expect_equal(unname(result["EuclideanD"]), sqrt(sum(colMeans(difference)^2)))
  expect_warning({ reordered = repeated_measures_test(a, b[n:1, ]) }, "used for matching")
  expect_equal(reordered, result)
  expect_warning({ positional = repeated_measures_test(a, b, rnames = FALSE) }, "same order")
  expect_equal(positional, result)
  expect_warning({ swapped = repeated_measures_test(b, a) }, "used for matching")
  expect_equal(swapped, result)
})

test_that("univariate repeated measures equal the paired t test", {
  blocks = statistical_blocks()
  a = blocks$x[, 1, drop = FALSE]
  b = blocks$y[, 1, drop = FALSE]
  expect_warning({ result = repeated_measures_test(a, b, rnames = FALSE) }, "same order")
  reference = t.test(a[, 1], b[, 1], paired = TRUE)
  expect_equal(unname(result["HotellingT2"]), unname(reference$statistic^2))
  expect_equal(unname(result["Fstat"]), unname(reference$statistic^2))
  expect_equal(unname(result["p_value"]), reference$p.value)
  expect_error(suppressWarnings(repeated_measures_test(a[-1, , drop = FALSE], b, rnames = FALSE)), "different number of rows")
  expect_error(suppressWarnings(repeated_measures_test(blocks$x, blocks$y[, 1:2], rnames = FALSE)), "different number of columns")
})

test_that("paired arrays agree with vectorized matrices and match specimen names", {
  skip_if_not_installed("Morpho")
  withr::local_seed(917)
  a = array(rnorm(96), dim = c(2, 2, 24),
            dimnames = list(NULL, NULL, paste0("specimen_", seq_len(24))))
  b = a + array(rnorm(96, mean = 0.3), dim = dim(a), dimnames = dimnames(a))
  matrix_a = Morpho::vecx(a)
  matrix_b = Morpho::vecx(b)
  expect_warning({ reference = repeated_measures_test(matrix_a, matrix_b) }, "used for matching")
  expect_warning({ actual = repeated_measures_test(a, b) }, "used for matching")
  expect_equal(actual, reference)
  expect_warning({ reordered = repeated_measures_test(a, b[, , 24:1]) }, "used for matching")
  expect_equal(reordered, reference)
  expect_warning({ positional = repeated_measures_test(a, b, rnames = FALSE) }, "same order")
  expect_equal(positional, reference)
  expect_warning({ mixed = repeated_measures_test(matrix_a, b) }, "used for matching")
  expect_equal(mixed, reference)
})

test_that("scaled eigenvalue variance agrees with a known covariance spectrum", {
  orthogonal = qr.Q(qr(contr.helmert(12)))[, 1:3]
  x = sweep(orthogonal, 2, sqrt(c(1, 4, 9) * 11), "*")
  result = scaled_variance_of_eigenvalues(x, boot = 0)
  expected = (3 * sum(c(1, 4, 9)^2) / sum(c(1, 4, 9))^2 - 1) / 2
  expect_equal(unname(result), c(expected, 3))
  isotropic = scaled_variance_of_eigenvalues(orthogonal, boot = 0)
  expect_equal(unname(isotropic), c(0, 3), tolerance = 1e-12)
  rotation = qr.Q(qr(matrix(c(1, 2, 3, 2, 5, 1, 4, 1, 2), 3)))
  transformed = sweep(7 * x %*% rotation, 2, c(-4, 8, 2), "+")
  expect_equal(scaled_variance_of_eigenvalues(as.data.frame(transformed), boot = 0), result)
  expect_error(scaled_variance_of_eigenvalues(list(1, 2), boot = 0), "matrix or a data frame")
})

test_that("eigenvalue bootstrap and rarefaction summaries describe their replicates", {
  withr::local_seed(831)
  x = statistical_blocks()$x
  result = scaled_variance_of_eigenvalues(x, boot = 40, rarefy = 12)
  expect_equal(dim(result$boot_sveig), c(40L, 2L))
  expect_equal(result$mean_boot_sveig, mean(result$boot_sveig[, 1]))
  expect_equal(result$median_dimensions_used, median(result$boot_sveig[, 2]))
  expect_length(result$boot_95ci_sveig, 4)
  expect_true(all(is.finite(result$boot_sveig)))
  expect_true(all(result$boot_sveig[, 1] >= 0 & result$boot_sveig[, 1] <= 1))
  expect_true(all(result$boot_sveig[, 2] <= ncol(x)))
  expect_true(result$boot_95ci_sveig[1] <= result$boot_95ci_sveig[2])
  expect_true(result$boot_95ci_sveig[3] <= result$boot_95ci_sveig[4])
  expect_equal(result$boot_95ci_sveig,
               c(sort(result$boot_sveig[, 1])[c(1, 39)], sort(result$boot_sveig[, 2])[c(1, 39)]))
  set.seed(831)
  expect_equal(scaled_variance_of_eigenvalues(x, boot = 40, rarefy = 12), result)
  expect_warning({ single = scaled_variance_of_eigenvalues(x, boot = 0, rarefy = 12) }, "SINGLE rarefaction")
  expect_true(single[1] >= 0 & single[1] <= 1)
  expect_equal(unname(single[2]), 3)
})

test_that("small bootstrap counts retain both confidence endpoints for both statistics", {
  withr::local_seed(918)
  x = statistical_blocks()$x
  for (boot in c(1L, 10L, 20L)) {
    result = scaled_variance_of_eigenvalues(x, boot = boot)
    expect_equal(dim(result$boot_sveig), c(boot, 2L))
    expect_length(result$boot_95ci_sveig, 4)
    expect_true(all(is.finite(result$boot_95ci_sveig)))
    # For these counts the rounded tail ranks, bounded to available draws, are the extremes.
    expect_equal(unname(result$boot_95ci_sveig),
                 c(range(result$boot_sveig[, 1]), range(result$boot_sveig[, 2])))
  }
})

test_that("shrinkage eigenvalue variance agrees with covariance concentration", {
  skip_if_not_installed("nlshrink")
  withr::local_seed(932)
  x = statistical_blocks()$x
  covariance = nlshrink::linshrink_cov(x)
  p = ncol(x)
  concentration = (p * sum(covariance^2) / sum(diag(covariance))^2 - 1) / (p - 1)
  result = scaled_variance_of_eigenvalues(x, boot = 0, shrinkage = TRUE)
  expect_equal(unname(result), c(concentration, p))
  bootstrap = scaled_variance_of_eigenvalues(x, boot = 40, shrinkage = TRUE)
  expect_equal(dim(bootstrap$boot_sveig), c(40L, 2L))
  expect_true(all(is.finite(bootstrap$boot_sveig)))
  expect_equal(bootstrap$mean_boot_sveig, mean(bootstrap$boot_sveig[, 1]))
  expect_warning({ rarefied = scaled_variance_of_eigenvalues(x, boot = 0, rarefy = 12, shrinkage = TRUE) }, "SINGLE rarefaction")
  expect_true(rarefied[1] >= 0 && rarefied[1] <= 1)
  expect_equal(unname(rarefied[2]), p)
})

test_that("Escoufier RV equals the cosine of centered specimen Gram matrices", {
  blocks = statistical_blocks()
  x = blocks$x
  y = blocks$y[, 1:2]
  gram_x = tcrossprod(scale(x, scale = FALSE))
  gram_y = tcrossprod(scale(y, scale = FALSE))
  reference = sum(gram_x * gram_y) / sqrt(sum(gram_x^2) * sum(gram_y^2))
  result = EscoufierRV(x, y)
  expect_equal(result, reference)
  expect_true(result >= 0 && result <= 1)
  expect_equal(EscoufierRV(y, x), result)
  expect_equal(EscoufierRV(x, x), 1)
  expect_equal(EscoufierRV(x * -3 + 5, y * 2 - 10), result)
  orthogonal = qr.Q(qr(contr.helmert(12)))
  expect_equal(EscoufierRV(orthogonal[, 1:2], orthogonal[, 3:4]), 0, tolerance = 1e-12)
  expect_error(EscoufierRV(x[-1, ], y), "same number of rows")
})

