# Test suite for RVcomparison function
# This test suite covers the major functionality of RVcomparison function
# including input validation, basic functionality, multi-group comparisons, 
# output structure, and reproducibility

# Test data preparation
setup_rvcomparison_test_data = function() {
  set.seed(123)
  
  # Create sample data for testing
  n_per_group = 15
  n_vars_block1 = 5
  n_vars_block2 = 4
  
  # Group A data
  Block1_A = matrix(rnorm(n_per_group * n_vars_block1, mean = 0), 
                    nrow = n_per_group, ncol = n_vars_block1)
  Block2_A = matrix(rnorm(n_per_group * n_vars_block2, mean = 0), 
                    nrow = n_per_group, ncol = n_vars_block2)
  
  # Group B data (slightly different mean)
  Block1_B = matrix(rnorm(n_per_group * n_vars_block1, mean = 0.5), 
                    nrow = n_per_group, ncol = n_vars_block1)
  Block2_B = matrix(rnorm(n_per_group * n_vars_block2, mean = 0.5), 
                    nrow = n_per_group, ncol = n_vars_block2)
  
  # Group C data 
  Block1_C = matrix(rnorm(n_per_group * n_vars_block1, mean = -0.3), 
                    nrow = n_per_group, ncol = n_vars_block1)
  Block2_C = matrix(rnorm(n_per_group * n_vars_block2, mean = -0.3), 
                    nrow = n_per_group, ncol = n_vars_block2)
  
  # Combine data for two-group tests
  Block1_2grp = rbind(Block1_A, Block1_B)
  Block2_2grp = rbind(Block2_A, Block2_B)
  groups_2grp = factor(c(rep("A", n_per_group), rep("B", n_per_group)))
  
  # Combine data for three-group tests
  Block1_3grp = rbind(Block1_A, Block1_B, Block1_C)
  Block2_3grp = rbind(Block2_A, Block2_B, Block2_C)
  groups_3grp = factor(c(rep("A", n_per_group), rep("B", n_per_group), rep("C", n_per_group)))
  
  list(
    Block1_2grp = Block1_2grp,
    Block2_2grp = Block2_2grp,
    groups_2grp = groups_2grp,
    Block1_3grp = Block1_3grp,
    Block2_3grp = Block2_3grp,
    groups_3grp = groups_3grp
  )
}

test_that("RVcomparison basic functionality works with two groups", {
  test_data = setup_rvcomparison_test_data()
  
  # Test basic two-group comparison
  result = RVcomparison(test_data$Block1_2grp, test_data$Block2_2grp, 
                       group = test_data$groups_2grp, perm = 99)
  
  # Should return a data frame
  expect_s3_class(result, "data.frame")
  
  # Should have exactly one row for one pairwise comparison
  expect_equal(nrow(result), 1)
  
  # Should have the correct groups
  expect_equal(result$group1, "A")
  expect_equal(result$group2, "B")
  
  # RV values should be numeric and between 0 and 1
  expect_true(is.numeric(result$Observed_RV_group1))
  expect_true(is.numeric(result$Observed_RV_group2))
  expect_true(result$Observed_RV_group1 >= 0 && result$Observed_RV_group1 <= 1)
  expect_true(result$Observed_RV_group2 >= 0 && result$Observed_RV_group2 <= 1)
  
  # P-value should be between 0 and 1
  expect_true(result$p_value >= 0 && result$p_value <= 1)
})

test_that("RVcomparison input validation works correctly", {
  test_data = setup_rvcomparison_test_data()
  
  # Test error for negative perm
  expect_error(RVcomparison(test_data$Block1_2grp, test_data$Block2_2grp, 
                           group = test_data$groups_2grp, perm = -5),
               "perm must be a positive integer")
  
  # Test error for non-integer perm
  expect_error(RVcomparison(test_data$Block1_2grp, test_data$Block2_2grp, 
                           group = test_data$groups_2grp, perm = 2.5),
               "perm must be a positive integer")
  
  # Test error for NA values in Block1
  Block1_with_NA = test_data$Block1_2grp
  Block1_with_NA[1, 1] = NA
  expect_error(RVcomparison(Block1_with_NA, test_data$Block2_2grp, 
                           group = test_data$groups_2grp, perm = 99),
               "Block1 contains NA values")
  
  # Test error for non-numeric data
  Block1_char = test_data$Block1_2grp
  Block1_char[, 1] = as.character(Block1_char[, 1])
  expect_error(RVcomparison(Block1_char, test_data$Block2_2grp, 
                           group = test_data$groups_2grp, perm = 99),
               "Block1 must contain only numeric variables")
  
  # Test error for single group
  groups_single = factor(rep("A", nrow(test_data$Block1_2grp)))
  expect_error(RVcomparison(test_data$Block1_2grp, test_data$Block2_2grp, 
                           group = groups_single, perm = 99),
               "At least two groups are required for RV comparison")
})

test_that("RVcomparison handles multiple groups with pairwise comparisons", {
  test_data = setup_rvcomparison_test_data()
  
  # Test three-group comparison
  result = RVcomparison(test_data$Block1_3grp, test_data$Block2_3grp, 
                       group = test_data$groups_3grp, perm = 99)
  
  # Should have 3 rows for 3 pairwise comparisons (A-B, A-C, B-C)
  expect_equal(nrow(result), 3)
  
  # Check that all group pairs are present
  group_pairs = paste(result$group1, result$group2, sep = "-")
  expected_pairs = c("A-B", "A-C", "B-C")
  expect_true(all(expected_pairs %in% group_pairs))
  
  # Each comparison should have valid RV values
  expect_true(all(result$Observed_RV_group1 >= 0 & result$Observed_RV_group1 <= 1))
  expect_true(all(result$Observed_RV_group2 >= 0 & result$Observed_RV_group2 <= 1))
  
  # All p-values should be valid
  expect_true(all(result$p_value >= 0 & result$p_value <= 1))
})

test_that("RVcomparison output has correct structure and column names", {
  test_data = setup_rvcomparison_test_data()
  
  result = RVcomparison(test_data$Block1_2grp, test_data$Block2_2grp, 
                       group = test_data$groups_2grp, perm = 99)
  
  # Check column names
  expected_cols = c("group1", "group2", "Observed_RV_group1", "Observed_RV_group2", 
                   "Absolute_difference_in_RV", "p_value")
  expect_equal(names(result), expected_cols)
  
  # Check data types
  expect_true(is.character(result$group1))
  expect_true(is.character(result$group2))
  expect_true(is.numeric(result$Observed_RV_group1))
  expect_true(is.numeric(result$Observed_RV_group2))
  expect_true(is.numeric(result$Absolute_difference_in_RV))
  expect_true(is.numeric(result$p_value))
  
  # Check that absolute difference is calculated correctly
  expected_diff = abs(result$Observed_RV_group1 - result$Observed_RV_group2)
  expect_equal(result$Absolute_difference_in_RV, expected_diff, tolerance = 1e-10)
})

test_that("RVcomparison produces reproducible results with same seed", {
  test_data = setup_rvcomparison_test_data()
  
  # Run the same test twice with same data and seed
  set.seed(456)
  result1 = RVcomparison(test_data$Block1_2grp, test_data$Block2_2grp, 
                        group = test_data$groups_2grp, perm = 99)
  
  set.seed(456)
  result2 = RVcomparison(test_data$Block1_2grp, test_data$Block2_2grp, 
                        group = test_data$groups_2grp, perm = 99)
  
  # Results should be identical
  expect_equal(result1, result2)
  
  # Observed RV values should be the same (they don't depend on permutations)
  expect_equal(result1$Observed_RV_group1, result2$Observed_RV_group1)
  expect_equal(result1$Observed_RV_group2, result2$Observed_RV_group2)
  
  # P-values should be the same with same seed
  expect_equal(result1$p_value, result2$p_value)
})