test_that("adjusted Rand tests recover known partition agreements", {
  withr::local_seed(601)
  a = c(1, 1, 2, 2)
  crossing = adjRand_test(a, c("a", "b", "a", "b"), perm = 39)
  matching = adjRand_test(rep(1:3, each = 3), rep(letters[1:3], each = 3),
                          perm = 39)

  expect_named(crossing, c("Adjusted_Rand_index", "p_value"))
  expect_equal(unname(crossing[1]), -0.5)
  expect_equal(unname(matching[1]), 1)
  expect_true(all(c(crossing[2], matching[2]) >= 1 / 40))
  expect_true(all(c(crossing[2], matching[2]) <= 1))
  expect_equal(unname(crossing[2] * 40), round(unname(crossing[2] * 40)))
  expect_error(adjRand_test(a, a[-1], perm = 39), "same length")
})

test_that("bootstrap disparity tails agree with returned replicates", {
  withr::local_seed(602)
  reference = matrix(c(0, 0, 2, 0, 0, 2, 2, 2), ncol = 2, byrow = TRUE)
  test = reference / 2
  result = BTailTest(reference, test, boot = 40)
  expected_labels = c("Multivariate variance", "Mean pairwise Euclidean distance")
  expect_named(result, c("BootstrappedSamplesEstimates", "pvalues"))
  expect_equal(dim(result$BootstrappedSamplesEstimates), c(40L, 2L))
  expect_equal(colnames(result$BootstrappedSamplesEstimates), expected_labels)
  expect_equal(rownames(result$pvalues), expected_labels)
  expect_equal(colnames(result$pvalues), c("Upper tail", "Lower tail", "Two-tail"))
  expect_true(all(is.finite(result$BootstrappedSamplesEstimates)))
  expect_true(all(result$BootstrappedSamplesEstimates >= 0))

  observed = c(2 / 3, (4 + 2 * sqrt(2)) / 6) # Unit square: four edges and two diagonals.
  for (i in seq_along(observed)) {
    replicates = result$BootstrappedSamplesEstimates[, i]
    tails = c(mean(replicates > observed[i]), mean(replicates < observed[i]))
    expect_equal(unname(result$pvalues[i, ]), c(tails, 2 * min(tails)))
  }
})

test_that("bootstrap disparity is invariant to separate group translations", {
  withr::local_seed(603)
  reference = matrix(c(0, 0, 2, 0, 0, 2, 2, 2), ncol = 2, byrow = TRUE)
  set.seed(603)
  original = BTailTest(reference, reference / 2, boot = 40)
  set.seed(603)
  translated = BTailTest(sweep(reference, 2, c(4, -8), "+"),
                         sweep(reference / 2, 2, c(-2, 6), "+"), boot = 40)
  expect_equal(translated, original)
})

test_that("critical angles match three-dimensional spherical cap geometry", {
  skip_if_not_installed("Morpho")
  alpha = c(0.01, 0.05, 0.1)
  angles = vapply(alpha, function(a) critical_angle(3, alpha = a), numeric(1))
  expect_equal(angles, acos(1 - 2 * alpha) * 180 / pi)
  expect_true(all(diff(angles) > 0))
})

test_that("critical angles validate their documented arguments", {
  expect_error(critical_angle(0), "positive integer")
  expect_error(critical_angle(2.5), "positive integer")
  expect_error(critical_angle("3"), "positive integer")
  expect_error(critical_angle(3, alpha = 0), "between 0 and 1")
  expect_error(critical_angle(3, alpha = 1), "between 0 and 1")
  expect_error(critical_angle(3, prec = "unknown"), "normal.*mpfr")
})

test_that("high precision critical angles agree at ordinary dimensions", {
  skip_if_not_installed("Morpho")
  skip_if_not_installed("Rmpfr")
  expect_equal(critical_angle(3, prec = "mpfr"), critical_angle(3), tolerance = 1e-10)
})

test_that("angle tests recover geometry and respect axis flipping", {
  skip_if_not_installed("Morpho")
  x = c(1, 0, 0)
  y = c(1 / 2, sqrt(3) / 2, 0)
  result = TestOfAngle(x, y)
  expect_true(all(c("angle", "angle_deg", "p.value", "critical_angle") %in% names(result)))
  expect_equal(result$angle, pi / 3)
  expect_equal(result$angle_deg, 60)
  expect_equal(result$p.value, 0.25)
  expect_equal(result$critical_angle, acos(0.9) * 180 / pi)
  expect_equal(TestOfAngle(2 * x, 3 * y), result)
  expect_equal(TestOfAngle(x, -y, flip = TRUE), result)
  expect_equal(TestOfAngle(x, -y, flip = FALSE)$angle_deg, 120)
  expect_equal(TestOfAngle(x, x)$angle_deg, 0)
})

test_that("disparity tests return known observed square statistics", {
  withr::local_seed(604)
  x = matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)
  result = disparity_test(x, 2 * x, perm = 39)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("Observed_grp1", "Observed_grp2", "difference", "p_value"))
  expect_equal(rownames(result), c("Multivariate variance", "Mean pairwise Euclidean distance"))
  expected = c(2 / 3, (4 + 2 * sqrt(2)) / 6)
  expect_equal(result$Observed_grp1, expected)
  expect_equal(result$Observed_grp2, expected * c(4, 2))
  expect_equal(result$difference, expected * c(3, 1))
  expect_true(all(result$p_value >= 1 / 40 & result$p_value <= 1))
  expect_equal(result$p_value * 40, round(result$p_value * 40))

  identical_groups = disparity_test(x, x, perm = 39)
  expect_equal(identical_groups$difference, c(0, 0))
  expect_equal(identical_groups$p_value, c(1, 1))
})

test_that("bootstrap mean distances recover an exact constant separation", {
  withr::local_seed(605)
  a = matrix(0, nrow = 4, ncol = 2)
  b = matrix(rep(c(3, 4), 6), ncol = 2, byrow = TRUE)
  result = dist_mean_boot(a, b, boot = 40, nA = 3, nB = 5)
  expect_named(result, c("Mean distance", "Median distance", "Lower CI extreme", "Upper CI extreme"))
  expect_equal(unname(result), rep(5, 4))
  expect_equal(dist_mean_boot(as.data.frame(a), as.data.frame(b), boot = 40), result)
})

test_that("bootstrap mean distances scale with coordinate units", {
  withr::local_seed(606)
  a = matrix(seq_len(16), ncol = 2)
  b = matrix(seq_len(12), ncol = 2)
  set.seed(606)
  original = dist_mean_boot(a, b, boot = 40)
  set.seed(606)
  scaled = dist_mean_boot(3 * a, 3 * b, boot = 40)
  expect_equal(scaled, 3 * original)
  expect_true(all(is.finite(original)))
  expect_true(original[3] <= original[2] && original[2] <= original[4])
})

test_that("bootstrap mean distances honor both requested resample sizes", {
  withr::local_seed(608)
  a = matrix(c(0, 0, 1, 3, 4, 2, 7, 1, 3, 8), ncol = 2, byrow = TRUE)
  b = matrix(c(2, 1, 6, 4, 3, 9, 8, 7), ncol = 2, byrow = TRUE)
  result = dist_mean_boot(a, b, boot = 80, ci = 0.8, nA = 3, nB = 7)

  set.seed(608)
  distances = vapply(seq_len(80), function(i) {
    weights_a = tabulate(sample.int(nrow(a), 3, replace = TRUE), nbins = nrow(a)) / 3
    weights_b = tabulate(sample.int(nrow(b), 7, replace = TRUE), nbins = nrow(b)) / 7
    difference = drop(crossprod(weights_a, a) - crossprod(weights_b, b))
    sqrt(sum(difference^2)) # Weighted centroids avoid duplicating row-extraction logic.
  }, numeric(1))
  expected = c(mean(distances), median(distances),
               quantile(distances, c(0.1, 0.9), type = 1, names = FALSE))
  expect_equal(unname(result), expected)

  set.seed(608)
  defaults = dist_mean_boot(a, b, boot = 80, ci = 0.8)
  expect_false(isTRUE(all.equal(unname(defaults), expected)))
})

test_that("bootstrap mean distances reject incompatible inputs", {
  x = matrix(seq_len(12), ncol = 2)
  expect_error(dist_mean_boot(x, x, ci = 0), "confidence interval")
  expect_error(dist_mean_boot(x, x, ci = 1), "confidence interval")
  expect_error(dist_mean_boot(x[1, , drop = FALSE], x), "more than one")
  expect_error(dist_mean_boot(x[, 1, drop = FALSE], x), "more than one")
  expect_error(dist_mean_boot(x, matrix(seq_len(18), ncol = 3)), "number of columns")
})

test_that("parallel analysis produces bounded permutation probabilities", {
  withr::local_seed(607)
  x = matrix(c(2, 4, 1, 8, 3, 9, 5, 7, 4, 6, 2, 1, 7, 3, 8, 5, 9, 6), ncol = 3)
  set.seed(607)
  result = parallel_analysis(x, perm = 39)
  expect_s3_class(result, "parallel_analysis")
  expect_length(result, 3)
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 1 / 40 & result <= 1))
  expect_equal(as.numeric(result) * 40, round(as.numeric(result) * 40))
  set.seed(607)
  expect_equal(parallel_analysis(as.data.frame(x), perm = 39), result)
  set.seed(607)
  expect_equal(parallel_analysis(x, perm = 39, fun = "fastSVD"), result)
  expect_equal(as.numeric(parallel_analysis(matrix(0, 6, 3), perm = 39)), rep(1, 3))
})

test_that("parallel analysis rejects inputs without multiple rows and columns", {
  expect_error(parallel_analysis(1:6, perm = 39), "matrix or a data.frame")
  expect_error(parallel_analysis(matrix(1:6, nrow = 1), perm = 39), "multiple rows and columns")
  expect_error(parallel_analysis(matrix(1:6, ncol = 1), perm = 39), "multiple rows and columns")
})

test_that("parallel analysis summary counts only initial significant axes", {
  cases = list(c(0.2, 0.3, 0.4), c(0.01, 0.05, 0.02), c(0.01, 0.2, 0.01), c(0.2, 0.01))
  expected = c(0, 3, 1, 0)
  for (i in seq_along(cases)) {
    object = structure(cases[[i]], class = c("parallel_analysis", "numeric"))
    output = capture.output(summary(object))
    expect_true(any(grepl(paste0("Number of tested axes: ", length(object)), output, fixed = TRUE)))
    expect_true(any(grepl(paste0("Number of axes significant at the 5% level: ", expected[i]), output, fixed = TRUE)))
  }
})

test_that("retired Kmultparallel entry points give defunct errors", {
  expect_error(Kmultparallel(), "geomorph::physignal", class = "defunctError")
  object = structure(list(), class = "parallel_Kmult")
  expect_error(print(object), "defunct", class = "defunctError")
  expect_error(plot(object), "defunct", class = "defunctError")
  expect_error(summary(object), "defunct", class = "defunctError")
})
