# Tests for GeometricMorphometricsMix Package

This directory contains comprehensive unit tests for the `disparity_resample()` function using the testthat framework.

## Test Structure

### `testthat.R`
Main test runner file that sets up the testing environment.

### Test Files

1. **`test-disparity_resample.R`** - Main comprehensive test suite
   - Input validation tests
   - Data type handling (matrix, data frame, vector, 3D array)
   - All statistical methods (multivariate variance, mean pairwise distance, etc.)
   - Bootstrap vs rarefaction methods
   - Confidence interval calculations
   - Missing data handling
   - Edge cases and error conditions

2. **`test-disparity_resample_helpers.R`** - Internal function tests
   - Tests for helper functions like `muvar()`, `meanpairwiseEuclideanD()`
   - Validation function tests
   - Data preparation function tests
   - Sample size resolution tests

3. **`test-disparity_resample_methods.R`** - S3 method tests
   - `print.disparity_resample()` method tests
   - `plot.disparity_resample()` method tests
   - Output formatting and display tests
   - Confidence interval overlap assessment tests

4. **`test-disparity_resample_integration.R`** - Integration and edge case tests
   - Tests with package datasets
   - Large dataset performance tests
   - Extreme value handling
   - Memory efficiency tests
   - Cross-platform compatibility tests

## Test Coverage

The test suite covers:

### Core Functionality
- ✅ All input data types (matrix, data.frame, vector, 3D array)
- ✅ All statistical methods
- ✅ Bootstrap and rarefaction resampling
- ✅ Single and multiple group analyses
- ✅ Custom sample sizes
- ✅ Confidence interval calculations

### Input Validation
- ✅ Parameter validation (CI, n_resamples, statistic, etc.)
- ✅ Data validation (missing data, empty datasets, etc.)
- ✅ Group validation (minimum sizes, matching lengths, etc.)
- ✅ Rarefaction-specific validation

### Statistical Accuracy
- ✅ Numerical correctness verification
- ✅ Reproducibility with set.seed()
- ✅ Edge cases (constant data, extreme values, etc.)

### S3 Methods
- ✅ Print method formatting and content
- ✅ Plot method ggplot2 integration
- ✅ Confidence interval overlap assessment
- ✅ Range vs confidence interval display

### Dependencies
- ✅ Optional package handling (geometry, nlshrink, Morpho, etc.)
- ✅ Graceful degradation when packages unavailable

## Running the Tests

To run all tests:
```r
library(testthat)
library(GeometricMorphometricsMix)
test_check("GeometricMorphometricsMix")
```

To run specific test files:
```r
test_file("tests/testthat/test-disparity_resample.R")
```

## Test Dependencies

The tests require:
- **testthat** (>= 3.0.0) - Core testing framework
- **GeometricMorphometricsMix** - The package being tested

Optional packages (tests will skip if not available):
- **MASS** - For multivariate normal data generation
- **geometry** - For convex hull volume tests
- **nlshrink** - For Claramunt proper variance tests
- **Morpho** - For 3D array handling tests
- **ggplot2** - For plot method tests

## Notes

- Tests use `skip_if_not_installed()` to gracefully handle missing optional packages
- Performance tests are marked with `skip_on_cran()` to avoid CI timeouts
- All tests use reproducible random seeds for consistency
- Error messages are tested to ensure they are informative

## Test Statistics

Total test count: ~100+ individual test cases across 4 test files
Coverage: All major functions and edge cases of `disparity_resample()`
