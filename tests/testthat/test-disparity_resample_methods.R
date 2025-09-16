# Test S3 methods for disparity_resample objects
# Tests for print.disparity_resample and plot.disparity_resample methods

setup_test_results = function() {
  set.seed(123)
  X1 = matrix(rnorm(20 * 3, mean = 0, sd = 1), nrow = 20, ncol = 3)
  X2 = matrix(rnorm(30 * 3, mean = 1, sd = 1), nrow = 30, ncol = 3)
  test_data = rbind(X1, X2)
  test_groups = factor(c(rep("A", 20), rep("B", 30)))
  
  # Create different types of results for testing
  multi_group_result = disparity_resample(test_data, group = test_groups, n_resamples = 50)
  single_group_result = disparity_resample(test_data, group = NULL, n_resamples = 50)
  range_result = disparity_resample(test_data, group = test_groups, n_resamples = 50, CI = 1)
  
  list(
    multi_group = multi_group_result,
    single_group = single_group_result,
    range = range_result
  )
}

test_that("print.disparity_resample works with multiple groups", {
  results = setup_test_results()
  
  # Should contain key components
  expect_output(print(results$multi_group), "Disparity resampling results")
  expect_output(print(results$multi_group), "Statistic: Multivariate variance")
  expect_output(print(results$multi_group), "Confidence level: 95%")
  expect_output(print(results$multi_group), "Confidence interval overlap assessment")
  
  # Should print the results table
  expect_output(print(results$multi_group), "group")
  expect_output(print(results$multi_group), "average")
  expect_output(print(results$multi_group), "CI_min")
  expect_output(print(results$multi_group), "CI_max")
})

test_that("print.disparity_resample works with single group", {
  results = setup_test_results()
  
  # Should contain key components but no overlap assessment
  expect_output(print(results$single_group), "Disparity resampling results")
  expect_output(print(results$single_group), "Statistic: Multivariate variance")
  expect_output(print(results$single_group), "Confidence level: 95%")
  
  # Should NOT contain overlap assessment for single group
  expect_output(print(results$single_group), "All", fixed = TRUE)
  output = capture.output(print(results$single_group))
  expect_false(any(grepl("overlap assessment", output)))
})

test_that("print.disparity_resample handles range results correctly", {
  results = setup_test_results()
  
  # Should show range instead of confidence interval
  expect_output(print(results$range), "Range: Full range")
  expect_output(print(results$range), "Range overlap assessment")
  
  # Should NOT show confidence level percentage
  output = capture.output(print(results$range))
  expect_false(any(grepl("Confidence level: [0-9]+%", output)))
})

test_that("print.disparity_resample returns object invisibly", {
  results = setup_test_results()
  
  # Should return the object invisibly
  returned = print(results$multi_group)
  expect_identical(returned, results$multi_group)
})

test_that("plot.disparity_resample works with multiple groups", {
  skip_if_not_installed("ggplot2")
  
  results = setup_test_results()
  
  plot_obj = plot(results$multi_group)
  
  # Should return a ggplot object
  expect_s3_class(plot_obj, "ggplot")

  # Check plot structure in a ggplot2-version-agnostic way
  built = ggplot2::ggplot_build(plot_obj)
  expect_s3_class(built, "ggplot_built")
  # At least one layer worth of data should exist
  expect_true(length(built$data) >= 1)
  # Labels should be present on the built plot
  expect_true(is.list(built$plot$labels))
  expect_true(length(built$plot$labels) >= 1)
})

test_that("plot.disparity_resample works with single group", {
  skip_if_not_installed("ggplot2")
  
  results = setup_test_results()
  
  plot_obj = plot(results$single_group)
  
  # Should return a ggplot object
  expect_s3_class(plot_obj, "ggplot")

  # Should handle single group case appropriately (check via ggplot_build)
  built = ggplot2::ggplot_build(plot_obj)
  expect_s3_class(built, "ggplot_built")
  expect_true(length(built$data) >= 1)
})

test_that("plot.disparity_resample handles range results", {
  skip_if_not_installed("ggplot2")
  
  results = setup_test_results()
  
  plot_obj = plot(results$range)
  
  # Should return a ggplot object
  expect_s3_class(plot_obj, "ggplot")
})

test_that("plot.disparity_resample accepts additional arguments", {
  skip_if_not_installed("ggplot2")
  
  results = setup_test_results()
  
  # Should accept additional arguments without error
  plot_obj = plot(results$multi_group, title = "Custom Title")
  expect_s3_class(plot_obj, "ggplot")
})

# Test the overlap assessment logic in print method
test_that("print.disparity_resample overlap assessment works correctly", {
  # Create results with known behavior for overlap assessment
  set.seed(42)
  # Use very different means and lower sample sizes to increase chance of non-overlap
  X1 = matrix(rnorm(10 * 2, mean = -5, sd = 0.5), nrow = 10, ncol = 2)  # Low mean, low variance
  X2 = matrix(rnorm(10 * 2, mean = 5, sd = 0.5), nrow = 10, ncol = 2)   # High mean, low variance
  test_data = rbind(X1, X2)
  test_groups = factor(c(rep("A", 10), rep("B", 10)))
  
  result = disparity_resample(test_data, group = test_groups, n_resamples = 100)
  
  # Capture output to check overlap assessment is present
  output = capture.output(print(result))
  overlap_line = output[grepl("overlap assessment", output)]
  expect_true(length(overlap_line) > 0)
  
  # Check that overlap assessment mentions either overlapping or non-overlapping
  has_overlap_info = any(grepl("overlap|Overlap", output))
  expect_true(has_overlap_info)
})

test_that("print.disparity_resample handles different statistics", {
  set.seed(123)
  test_data = matrix(rnorm(50 * 3), nrow = 50, ncol = 3)
  test_groups = factor(c(rep("A", 25), rep("B", 25)))
  
  # Test different statistics
  result_mpd = disparity_resample(test_data, group = test_groups, 
                                  statistic = "mean_pairwise_euclidean_distance", 
                                  n_resamples = 20)
  
  expect_output(print(result_mpd), "Mean pairwise Euclidean distance")
})

test_that("print.disparity_resample handles univariate data", {
  set.seed(123)
  univar_data = rnorm(50)
  univar_groups = factor(c(rep("X", 25), rep("Y", 25)))
  
  result = disparity_resample(univar_data, group = univar_groups, n_resamples = 20)
  
  expect_output(print(result), "Statistic: Variance")
})

test_that("S3 method dispatch works correctly", {
  results = setup_test_results()
  
  # Check that the object has the correct class
  expect_s3_class(results$multi_group, "disparity_resample")
  expect_s3_class(results$multi_group, "list")
  
  # Check that S3 methods are properly dispatched
  expect_true(exists("print.disparity_resample"))
  expect_true(exists("plot.disparity_resample"))
})
