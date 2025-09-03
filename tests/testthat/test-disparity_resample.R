# Test suite for disparity_resample function
# This comprehensive test suite covers all major aspects of the disparity_resample function
# including input validation, different data types, statistics, and resampling methods

# Test data preparation
setup_test_data = function() {
  set.seed(123)
  
  # Basic multivariate data
  X1 = matrix(rnorm(20 * 5, mean = 0, sd = 1), nrow = 20, ncol = 5)
  X2 = matrix(rnorm(30 * 5, mean = 1, sd = 1), nrow = 30, ncol = 5)
  basic_data = rbind(X1, X2)
  basic_groups = factor(c(rep("A", 20), rep("B", 30)))
  
  # Univariate data
  univar_data = rnorm(50)
  univar_groups = factor(c(rep("X", 25), rep("Y", 25)))
  
  # Small dataset for testing edge cases
  small_data = matrix(rnorm(6 * 3), nrow = 6, ncol = 3)
  small_groups = factor(c(rep("G1", 3), rep("G2", 3)))
  
  # Data with missing values
  missing_data = basic_data
  missing_data[1, 1] = NA
  missing_data[25, 3] = NA
  
  list(
    basic_data = basic_data,
    basic_groups = basic_groups,
    univar_data = univar_data,
    univar_groups = univar_groups,
    small_data = small_data,
    small_groups = small_groups,
    missing_data = missing_data
  )
}

test_that("disparity_resample returns correct structure", {
  test_data = setup_test_data()
  
  result = disparity_resample(test_data$basic_data, group = test_data$basic_groups, 
                              n_resamples = 10)
  
  # Check return structure
  expect_s3_class(result, "disparity_resample")
  expect_type(result, "list")
  expect_named(result, c("chosen_statistic", "results", "resampled_values", "CI_level"))
  
  # Check results data frame structure
  expect_s3_class(result$results, "data.frame")
  expect_named(result$results, c("group", "average", "CI_min", "CI_max"))
  expect_equal(nrow(result$results), 2)  # Two groups
  
  # Check resampled values structure for multiple groups
  expect_type(result$resampled_values, "list")
  expect_length(result$resampled_values, 2)
  expect_named(result$resampled_values, c("A", "B"))
  expect_length(result$resampled_values$A, 10)
  expect_length(result$resampled_values$B, 10)
})

test_that("disparity_resample works with single group", {
  test_data = setup_test_data()
  
  result = disparity_resample(test_data$basic_data, group = NULL, n_resamples = 10)
  
  # Check single group structure
  expect_equal(nrow(result$results), 1)
  expect_equal(result$results$group, "All")
  
  # Check resampled values structure for single group
  expect_type(result$resampled_values, "double")
  expect_length(result$resampled_values, 10)
})

test_that("disparity_resample works with univariate data", {
  test_data = setup_test_data()
  
  result = disparity_resample(test_data$univar_data, group = test_data$univar_groups, 
                              n_resamples = 10)
  
  expect_equal(result$chosen_statistic, "Variance")
  expect_equal(nrow(result$results), 2)
  expect_true(all(result$results$average > 0))  # Variances should be positive
})

test_that("disparity_resample validates input parameters correctly", {
  test_data = setup_test_data()
  
  # Test invalid CI values
  expect_error(disparity_resample(test_data$basic_data, CI = 0))
  expect_error(disparity_resample(test_data$basic_data, CI = 1.5))
  expect_error(disparity_resample(test_data$basic_data, CI = c(0.8, 0.9)))
  
  # Test invalid n_resamples
  expect_error(disparity_resample(test_data$basic_data, n_resamples = 0))
  expect_error(disparity_resample(test_data$basic_data, n_resamples = -5))
  expect_error(disparity_resample(test_data$basic_data, n_resamples = 1.5))
  
  # Test invalid bootstrap_rarefaction
  expect_error(disparity_resample(test_data$basic_data, bootstrap_rarefaction = "invalid"))
  
  # Test invalid statistic
  expect_error(disparity_resample(test_data$basic_data, statistic = "invalid_stat"))
})

test_that("disparity_resample validates rarefaction requirements", {
  test_data = setup_test_data()
  
  # Rarefaction requires sample_size
  expect_error(disparity_resample(test_data$basic_data, group = test_data$basic_groups,
                                 bootstrap_rarefaction = "rarefaction"))
  
  # sample_size must be valid for rarefaction
  expect_error(disparity_resample(test_data$basic_data, group = test_data$basic_groups,
                                 bootstrap_rarefaction = "rarefaction", sample_size = 0))
  expect_error(disparity_resample(test_data$basic_data, group = test_data$basic_groups,
                                 bootstrap_rarefaction = "rarefaction", sample_size = 1.5))
  
  # sample_size cannot exceed smallest group size
  expect_error(disparity_resample(test_data$basic_data, group = test_data$basic_groups,
                                 bootstrap_rarefaction = "rarefaction", sample_size = 25))
})

test_that("disparity_resample validates group requirements", {
  test_data = setup_test_data()
  
  # Groups must have minimum size for variance calculation
  tiny_data = matrix(rnorm(2 * 3), nrow = 2, ncol = 3)
  tiny_groups = factor(c("A", "B"))
  expect_error(disparity_resample(tiny_data, group = tiny_groups))
  
  # Group length must match data
  wrong_groups = factor(c(rep("A", 10), rep("B", 15)))  # Wrong length
  expect_error(disparity_resample(test_data$basic_data, group = wrong_groups))
})

test_that("disparity_resample handles missing data correctly", {
  test_data = setup_test_data()
  
  expect_warning(result <- disparity_resample(test_data$missing_data, 
                                             group = test_data$basic_groups, 
                                             n_resamples = 10),
                 regexp = "missing|NA|removed")
  
  # Should still return valid results after removing missing data
  expect_s3_class(result, "disparity_resample")
  expect_equal(nrow(result$results), 2)
})

test_that("disparity_resample statistics work correctly", {
  test_data = setup_test_data()
  
  # Test multivariate variance
  result_mv = disparity_resample(test_data$basic_data, group = test_data$basic_groups,
                                 statistic = "multivariate_variance", n_resamples = 10)
  expect_equal(result_mv$chosen_statistic, "Multivariate variance")
  
  # Test mean pairwise Euclidean distance
  result_mpd = disparity_resample(test_data$basic_data, group = test_data$basic_groups,
                                  statistic = "mean_pairwise_euclidean_distance", n_resamples = 10)
  expect_equal(result_mpd$chosen_statistic, "Mean pairwise Euclidean distance")
})

test_that("disparity_resample convex hull volume works with adequate data", {
  skip_if_not_installed("geometry")
  
  # Create larger dataset for convex hull (needs n > p)
  set.seed(123)
  large_data = matrix(rnorm(100 * 3), nrow = 100, ncol = 3)  # 100 obs, 3 vars
  large_groups = factor(c(rep("A", 50), rep("B", 50)))
  
  result_cv = disparity_resample(large_data, group = large_groups,
                                 statistic = "convex_hull_volume", n_resamples = 10)
  expect_equal(result_cv$chosen_statistic, "Convex hull volume")
  expect_true(all(result_cv$results$average > 0))  # Volumes should be positive
})

test_that("disparity_resample convex hull volume validates dimensions", {
  skip_if_not_installed("geometry")
  
  # Create test data where convex hull would fail: few observations in high dimensions
  small_high_dim_data = matrix(rnorm(6 * 10), nrow = 6, ncol = 10)  # 6 obs, 10 dims
  small_groups = factor(c(rep("A", 3), rep("B", 3)))
  
  # Should fail when n <= p (insufficient observations for convex hull)
  expect_error(disparity_resample(small_high_dim_data, group = small_groups,
                                 statistic = "convex_hull_volume"),
               regexp = "convex hull|dimensions|insufficient|observations")
})

test_that("disparity_resample Claramunt proper variance works", {
  skip_if_not_installed("nlshrink")
  
  test_data = setup_test_data()
  
  result_cpv = disparity_resample(test_data$basic_data, group = test_data$basic_groups,
                                  statistic = "claramunt_proper_variance", n_resamples = 10)
  expect_equal(result_cpv$chosen_statistic, "Claramunt proper variance")
  expect_true(all(result_cpv$results$average > 0))  # Should be positive
})

test_that("disparity_resample bootstrap vs rarefaction work correctly", {
  test_data = setup_test_data()
  
  # Bootstrap (default)
  result_boot = disparity_resample(test_data$basic_data, group = test_data$basic_groups,
                                   bootstrap_rarefaction = "bootstrap", n_resamples = 10)
  
  # Rarefaction to smallest group
  result_rare = disparity_resample(test_data$basic_data, group = test_data$basic_groups,
                                   bootstrap_rarefaction = "rarefaction", 
                                   sample_size = "smallest", n_resamples = 10)
  
  expect_s3_class(result_boot, "disparity_resample")
  expect_s3_class(result_rare, "disparity_resample")
  
  # Both should have same structure
  expect_equal(names(result_boot), names(result_rare))
  expect_equal(nrow(result_boot$results), nrow(result_rare$results))
})

test_that("disparity_resample custom sample sizes work", {
  test_data = setup_test_data()
  
  # Bootstrap with custom sample size
  result_boot_custom = disparity_resample(test_data$basic_data, group = test_data$basic_groups,
                                          bootstrap_rarefaction = "bootstrap", 
                                          sample_size = 15, n_resamples = 10)
  
  # Rarefaction with custom sample size
  result_rare_custom = disparity_resample(test_data$basic_data, group = test_data$basic_groups,
                                          bootstrap_rarefaction = "rarefaction", 
                                          sample_size = 15, n_resamples = 10)
  
  expect_s3_class(result_boot_custom, "disparity_resample")
  expect_s3_class(result_rare_custom, "disparity_resample")
})

test_that("disparity_resample CI levels work correctly", {
  test_data = setup_test_data()
  
  # Test different CI levels
  result_95 = disparity_resample(test_data$basic_data, group = test_data$basic_groups,
                                 CI = 0.95, n_resamples = 100)
  result_99 = disparity_resample(test_data$basic_data, group = test_data$basic_groups,
                                 CI = 0.99, n_resamples = 100)
  result_full = disparity_resample(test_data$basic_data, group = test_data$basic_groups,
                                   CI = 1, n_resamples = 100)
  
  expect_equal(result_95$CI_level, 0.95)
  expect_equal(result_99$CI_level, 0.99)
  expect_equal(result_full$CI_level, 1)
  
  # 99% CI should be wider than 95% CI
  expect_true(all((result_99$results$CI_max - result_99$results$CI_min) >= 
                  (result_95$results$CI_max - result_95$results$CI_min)))
})

test_that("disparity_resample handles data frame input", {
  test_data = setup_test_data()
  
  # Convert to data frame
  df_data = as.data.frame(test_data$basic_data)
  
  result = disparity_resample(df_data, group = test_data$basic_groups, n_resamples = 10)
  expect_s3_class(result, "disparity_resample")
})

test_that("disparity_resample handles 3D array input", {
  skip_if_not_installed("Morpho")
  
  # Create 3D array (landmarks x dimensions x specimens)
  set.seed(123)
  landmarks_3d = array(rnorm(5 * 2 * 20), dim = c(5, 2, 20))  # 5 landmarks, 2D, 20 specimens
  groups_3d = factor(c(rep("A", 10), rep("B", 10)))
  
  result = disparity_resample(landmarks_3d, group = groups_3d, n_resamples = 10)
  expect_s3_class(result, "disparity_resample")
})

test_that("disparity_resample sample_size='smallest' validation works", {
  test_data = setup_test_data()
  
  # Should work with multiple groups
  result = disparity_resample(test_data$basic_data, group = test_data$basic_groups,
                              bootstrap_rarefaction = "rarefaction", 
                              sample_size = "smallest", n_resamples = 10)
  expect_s3_class(result, "disparity_resample")
  
  # Should fail with single group
  expect_error(disparity_resample(test_data$basic_data, group = NULL,
                                 bootstrap_rarefaction = "rarefaction", 
                                 sample_size = "smallest"))
})

test_that("print.disparity_resample method works", {
  test_data = setup_test_data()
  
  result = disparity_resample(test_data$basic_data, group = test_data$basic_groups, 
                              n_resamples = 10)
  
  # Should not throw error when printing
  expect_output(print(result), "Disparity resampling results")
  expect_output(print(result), "Multivariate variance")
  expect_output(print(result), "Confidence interval overlap assessment")
})

test_that("plot.disparity_resample method works", {
  skip_if_not_installed("ggplot2")
  
  test_data = setup_test_data()
  
  result = disparity_resample(test_data$basic_data, group = test_data$basic_groups, 
                              n_resamples = 10)
  
  # Should return a ggplot object
  p = plot(result)
  expect_s3_class(p, "ggplot")
})

test_that("disparity_resample numerical accuracy", {
  # Test with known data to verify statistical correctness
  set.seed(42)
  
  # Create data with known variance
  known_data = matrix(c(rep(0, 50), rep(1, 50)), nrow = 50, ncol = 2)
  var_x1 = var(known_data[, 1])  # Should be 0.25
  var_x2 = var(known_data[, 2])  # Should be 0.25
  expected_muvar = var_x1 + var_x2  # Should be 0.5
  
  result = disparity_resample(known_data, n_resamples = 1000)
  
  # Bootstrap average should be close to theoretical value
  expect_equal(result$results$average, expected_muvar, tolerance = 0.1)
})

test_that("disparity_resample reproducibility with set.seed", {
  test_data = setup_test_data()
  
  # Run with same seed twice
  set.seed(456)
  result1 = disparity_resample(test_data$basic_data, group = test_data$basic_groups, 
                               n_resamples = 10)
  
  set.seed(456)
  result2 = disparity_resample(test_data$basic_data, group = test_data$basic_groups, 
                               n_resamples = 10)
  
  # Results should be identical
  expect_equal(result1$resampled_values, result2$resampled_values)
  expect_equal(result1$results, result2$results)
})

test_that("disparity_resample edge case: minimum valid group sizes", {
  # Test with exactly minimum group size (2 observations per group)
  min_data = matrix(rnorm(4 * 3), nrow = 4, ncol = 3)
  min_groups = factor(c("A", "A", "B", "B"))
  
  result = disparity_resample(min_data, group = min_groups, n_resamples = 10)
  expect_s3_class(result, "disparity_resample")
  expect_equal(nrow(result$results), 2)
})
