# Test suite for EscoufierRVrarefy function
# This test suite covers the major functionality of RVrarefied function
# including input validation, single-group vs multi-group analysis, and S3 methods

# Test data preparation
setup_rv_test_data = function() {
  set.seed(42)
  
  # Basic test data - two blocks of variables
  n_obs = 30
  n_vars_block1 = 5
  n_vars_block2 = 4
  
  Block1 = matrix(rnorm(n_obs * n_vars_block1), nrow = n_obs, ncol = n_vars_block1)
  Block2 = matrix(rnorm(n_obs * n_vars_block2), nrow = n_obs, ncol = n_vars_block2)
  
  # Multi-group data
  n_total = 60
  Block1_multi = matrix(rnorm(n_total * n_vars_block1), nrow = n_total, ncol = n_vars_block1)
  Block2_multi = matrix(rnorm(n_total * n_vars_block2), nrow = n_total, ncol = n_vars_block2)
  groups = factor(c(rep("GroupA", 30), rep("GroupB", 30)))
  
  list(
    Block1 = Block1,
    Block2 = Block2,
    Block1_multi = Block1_multi,
    Block2_multi = Block2_multi,
    groups = groups
  )
}

test_that("RVrarefied returns correct structure for single-group analysis", {
  test_data = setup_rv_test_data()
  
  result = RVrarefied(test_data$Block1, test_data$Block2, reps = 10, samplesize = 20)
  
  # Check return structure
  expect_s3_class(result, "EscoufierRVrarefy")
  expect_type(result, "list")
  expect_named(result, c("results", "AllRarefiedSamples"))
  
  # Check results data frame structure
  expect_s3_class(result$results, "data.frame")
  expect_equal(nrow(result$results), 1)
  expect_named(result$results, c("group", "Mean", "Median", "CI_min", "CI_max"))
  expect_equal(result$results$group, "All")
  
  # Check AllRarefiedSamples structure
  expect_type(result$AllRarefiedSamples, "list")
  expect_named(result$AllRarefiedSamples, "All")
  expect_equal(length(result$AllRarefiedSamples$All), 10)
  
  # Check that values are numeric and reasonable
  expect_true(all(is.numeric(result$results$Mean)))
  expect_true(all(result$results$CI_min <= result$results$CI_max))
  expect_true(all(result$results$Mean >= 0))  # RV should be non-negative
})

test_that("RVrarefied returns correct structure for multi-group analysis", {
  test_data = setup_rv_test_data()
  
  result = RVrarefied(test_data$Block1_multi, test_data$Block2_multi, 
                     reps = 10, samplesize = 20, group = test_data$groups)
  
  # Check return structure
  expect_s3_class(result, "EscoufierRVrarefy")
  expect_type(result, "list")
  expect_named(result, c("results", "AllRarefiedSamples"))
  
  # Check results data frame structure for multiple groups
  expect_s3_class(result$results, "data.frame")
  expect_equal(nrow(result$results), 2)
  expect_named(result$results, c("group", "Mean", "Median", "CI_min", "CI_max"))
  expect_setequal(result$results$group, c("GroupA", "GroupB"))
  
  # Check AllRarefiedSamples structure for multiple groups
  expect_type(result$AllRarefiedSamples, "list")
  expect_named(result$AllRarefiedSamples, c("GroupA", "GroupB"))
  expect_equal(length(result$AllRarefiedSamples$GroupA), 10)
  expect_equal(length(result$AllRarefiedSamples$GroupB), 10)
})

test_that("RVrarefied input validation works correctly", {
  test_data = setup_rv_test_data()
  
  # Test mismatched Block dimensions
  Block1_wrong = matrix(rnorm(20 * 5), nrow = 20, ncol = 5)  # 20 rows instead of 30
  expect_error(RVrarefied(Block1_wrong, test_data$Block2, reps = 10, samplesize = 15),
               "same number of rows")
  
  # Test invalid CI values
  expect_error(RVrarefied(test_data$Block1, test_data$Block2, reps = 10, samplesize = 15, CI = 0),
               "CI must be between 0 and 1")
  expect_error(RVrarefied(test_data$Block1, test_data$Block2, reps = 10, samplesize = 15, CI = 1.5),
               "CI must be between 0 and 1")
  
  # Test invalid group length
  wrong_groups = factor(rep("A", 25))  # 25 instead of 30
  expect_error(RVrarefied(test_data$Block1, test_data$Block2, reps = 10, samplesize = 15, 
                         group = wrong_groups),
               "same length as the number of observations")
  
  # Test sample size too large
  expect_error(RVrarefied(test_data$Block1_multi, test_data$Block2_multi, reps = 10, 
                         samplesize = 35, group = test_data$groups),
               "fewer observations than the requested sample size")
})

test_that("RVrarefied CI parameter works correctly", {
  test_data = setup_rv_test_data()
  
  # Test different CI levels
  result_90 = RVrarefied(test_data$Block1_multi, test_data$Block2_multi, 
                        reps = 20, samplesize = 20, group = test_data$groups, CI = 0.90)
  result_99 = RVrarefied(test_data$Block1_multi, test_data$Block2_multi, 
                        reps = 20, samplesize = 20, group = test_data$groups, CI = 0.99)
  
  # 99% CI should be wider than 90% CI (generally)
  ci_width_90 = result_90$results$CI_max - result_90$results$CI_min
  ci_width_99 = result_99$results$CI_max - result_99$results$CI_min
  
  # At least one group should have wider CI for 99% (this is probabilistic)
  expect_true(any(ci_width_99 >= ci_width_90))
})

test_that("S3 methods work correctly", {
  test_data = setup_rv_test_data()
  
  result = RVrarefied(test_data$Block1_multi, test_data$Block2_multi, 
                     reps = 10, samplesize = 20, group = test_data$groups)
  
  # Test print method (should not error)
  expect_output(print(result), "Rarefied Escoufier RV coefficient")
  expect_output(print(result), "GroupA")
  expect_output(print(result), "GroupB")
  
  # Test plot method exists and returns ggplot (if ggplot2 is available)
  skip_if_not_installed("ggplot2")
  plot_result = plot(result)
  expect_s3_class(plot_result, "ggplot")
})
