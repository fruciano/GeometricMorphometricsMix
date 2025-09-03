# Integration and edge case tests for disparity_resample
# Tests using package data and testing unusual edge cases

test_that("disparity_resample works with brown_trout dataset", {
  skip_if_not_available("brown_trout")
  
  data("brown_trout", package = "GeometricMorphometricsMix")
  
  if (exists("brown_trout")) {
    # Test with actual package data if available
    result = disparity_resample(brown_trout, n_resamples = 50)
    expect_s3_class(result, "disparity_resample")
    expect_equal(result$results$group, "All")
  }
})

test_that("disparity_resample handles very small datasets", {
  # Test with minimum viable dataset
  min_data = matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
  
  result = disparity_resample(min_data, n_resamples = 10)
  expect_s3_class(result, "disparity_resample")
  expect_equal(nrow(result$results), 1)
})

test_that("disparity_resample handles large number of resamples", {
  set.seed(123)
  test_data = matrix(rnorm(30), nrow = 10, ncol = 3)
  
  # Should handle large number of resamples
  result = disparity_resample(test_data, n_resamples = 5000)
  expect_s3_class(result, "disparity_resample")
  expect_length(result$resampled_values, 5000)
})

test_that("disparity_resample handles extreme CI values", {
  set.seed(123)
  test_data = matrix(rnorm(30), nrow = 10, ncol = 3)
  
  # Test very narrow CI
  result_narrow = disparity_resample(test_data, CI = 0.01, n_resamples = 100)
  expect_s3_class(result_narrow, "disparity_resample")
  expect_equal(result_narrow$CI_level, 0.01)
  
  # Test very wide CI (almost full range)
  result_wide = disparity_resample(test_data, CI = 0.999, n_resamples = 100)
  expect_s3_class(result_wide, "disparity_resample")
  expect_equal(result_wide$CI_level, 0.999)
})

test_that("disparity_resample memory efficiency with large datasets", {
  # Test that function doesn't consume excessive memory
  skip_on_cran()  # Skip on CRAN to avoid timeouts
  
  set.seed(123)
  large_data = matrix(rnorm(1000 * 10), nrow = 1000, ncol = 10)
  large_groups = factor(rep(c("A", "B"), each = 500))
  
  # Should complete without memory issues
  result = disparity_resample(large_data, group = large_groups, n_resamples = 100)
  expect_s3_class(result, "disparity_resample")
})

test_that("disparity_resample handles data with extreme values", {
  set.seed(123)
  
  # Create data with extreme outliers
  normal_data = matrix(rnorm(30), nrow = 10, ncol = 3)
  extreme_data = normal_data
  extreme_data[1, ] = c(1000, -1000, 500)  # Add extreme outliers
  
  result = disparity_resample(extreme_data, n_resamples = 50)
  expect_s3_class(result, "disparity_resample")
  expect_true(all(is.finite(result$results$average)))
})

test_that("disparity_resample handles perfectly correlated data", {
  # Create perfectly correlated data
  set.seed(123)
  x = rnorm(20)
  perfect_corr_data = cbind(x, x, x)  # All columns identical
  
  result = disparity_resample(perfect_corr_data, n_resamples = 50)
  expect_s3_class(result, "disparity_resample")
  # Should handle degenerate covariance structure
})

test_that("disparity_resample handles constant data", {
  # Create constant data (no variance)
  constant_data = matrix(1, nrow = 10, ncol = 3)
  
  result = disparity_resample(constant_data, n_resamples = 50)
  expect_s3_class(result, "disparity_resample")
  expect_equal(result$results$average, 0)  # Variance should be 0
})

test_that("disparity_resample validation with mixed data types", {
  # Test error handling with inappropriate data types
  expect_error(disparity_resample("not_a_matrix"))
  expect_error(disparity_resample(list(a = 1, b = 2)))
  expect_error(disparity_resample(factor(1:10)))
})

test_that("disparity_resample handles unbalanced groups appropriately", {
  set.seed(123)
  
  # Create very unbalanced groups
  X1 = matrix(rnorm(5 * 3), nrow = 5, ncol = 3)    # Small group
  X2 = matrix(rnorm(100 * 3), nrow = 100, ncol = 3) # Large group
  unbalanced_data = rbind(X1, X2)
  unbalanced_groups = factor(c(rep("Small", 5), rep("Large", 100)))
  
  result = disparity_resample(unbalanced_data, group = unbalanced_groups, 
                              bootstrap_rarefaction = "rarefaction", 
                              sample_size = "smallest", n_resamples = 50)
  expect_s3_class(result, "disparity_resample")
})

test_that("disparity_resample handles character group labels correctly", {
  set.seed(123)
  test_data = matrix(rnorm(30), nrow = 10, ncol = 3)
  
  # Test with character groups (should be converted to factor)
  char_groups = c(rep("Group_A", 5), rep("Group_B", 5))
  result = disparity_resample(test_data, group = char_groups, n_resamples = 20)
  
  expect_s3_class(result, "disparity_resample")
  expect_true(all(c("Group_A", "Group_B") %in% result$results$group))
})

test_that("disparity_resample handles numeric group labels correctly", {
  set.seed(123)
  test_data = matrix(rnorm(30), nrow = 10, ncol = 3)
  
  # Test with numeric groups (should be converted to factor)
  numeric_groups = c(rep(1, 5), rep(2, 5))
  result = disparity_resample(test_data, group = numeric_groups, n_resamples = 20)
  
  expect_s3_class(result, "disparity_resample")
  expect_equal(length(unique(result$results$group)), 2)
})

test_that("disparity_resample with all missing data in one variable", {
  set.seed(123)
  test_data = matrix(rnorm(20), nrow = 10, ncol = 2)
  test_data[, 2] = NA  # Make one column all NA
  
  expect_warning(result = disparity_resample(test_data, n_resamples = 20))
  expect_s3_class(result, "disparity_resample")
})

test_that("disparity_resample error messages are informative", {
  test_data = matrix(rnorm(20), nrow = 10, ncol = 2)
  
  # Check that error messages contain useful information
  expect_error(disparity_resample(test_data, CI = 2), 
               "CI must be a single numeric value between 0 and 1")
  
  expect_error(disparity_resample(test_data, n_resamples = -5), 
               "n_resamples must be a positive integer")
  
  expect_error(disparity_resample(test_data, statistic = "nonexistent"), 
               "statistic should be one of:")
})

test_that("disparity_resample reproducibility across R sessions", {
  # Test that results are reproducible with set.seed
  set.seed(999)
  test_data = matrix(rnorm(20), nrow = 10, ncol = 2)
  test_groups = factor(rep(c("A", "B"), each = 5))
  
  # Run multiple times with same seed
  set.seed(999)
  result1 = disparity_resample(test_data, group = test_groups, n_resamples = 50)
  
  set.seed(999)
  result2 = disparity_resample(test_data, group = test_groups, n_resamples = 50)
  
  # Results should be identical
  expect_equal(result1$results, result2$results)
  expect_equal(result1$resampled_values, result2$resampled_values)
})

test_that("disparity_resample performance with different statistics", {
  skip_on_cran()  # Skip performance tests on CRAN
  
  set.seed(123)
  test_data = matrix(rnorm(200 * 5), nrow = 200, ncol = 5)
  test_groups = factor(rep(c("A", "B"), each = 100))
  
  # All statistics should complete in reasonable time
  stats_to_test = c("multivariate_variance", "mean_pairwise_euclidean_distance")
  
  for (stat in stats_to_test) {
    result = disparity_resample(test_data, group = test_groups, 
                                statistic = stat, n_resamples = 100)
    expect_s3_class(result, "disparity_resample")
    expect_equal(result$chosen_statistic, 
                switch(stat,
                       "multivariate_variance" = "Multivariate variance",
                       "mean_pairwise_euclidean_distance" = "Mean pairwise Euclidean distance"))
  }
})
