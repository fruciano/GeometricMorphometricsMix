# Suggested Improvements for disparity_resample()

## 1. MODULARIZATION (High Priority)

### A. Extract Input Validation Functions
```r
validate_disparity_inputs <- function(Data, group, n_resamples, statistic, CI, 
                                    bootstrap_rarefaction, sample_size) {
  # Move all validation logic here
  # Return standardized inputs or stop with clear errors
}

prepare_data_and_groups <- function(Data, group) {
  # Handle data coercion, missing data, group factor creation
  # Return list with cleaned data and group information
}

validate_group_sizes <- function(group_factor, statistic, sample_size_num = NULL) {
  # Check minimum group sizes for different statistics
  # Separate convex hull validation
}
```

### B. Extract Statistical Computation Functions
```r
create_statistic_registry <- function() {
  # Registry pattern for extensible statistics
  list(
    multivariate_variance = list(
      name = "Multivariate variance",
      fn = muvar,
      requires_packages = NULL,
      min_observations = 2
    ),
    mean_pairwise_euclidean_distance = list(
      name = "Mean pairwise Euclidean distance", 
      fn = meanpairwiseEuclideanD,
      requires_packages = NULL,
      min_observations = 2
    ),
    convex_hull_volume = list(
      name = "Convex hull volume",
      fn = function(X) geometry::convhulln(X, options="FA")$vol,
      requires_packages = "geometry",
      min_observations_formula = function(n_vars) n_vars + 1
    ),
    claramunt_proper_variance = list(
      name = "Claramunt proper variance",
      fn = claramunt_proper_variance,
      requires_packages = "nlshrink", 
      min_observations = 2
    )
  )
}

get_statistic_function <- function(statistic, original_input_is_vector) {
  # Return appropriate function based on statistic choice
}
```

### C. Extract Resampling Logic
```r
perform_bootstrap_resampling <- function(data_group, n_resamples, stat_fn) {
  # Bootstrap logic separated from main function
}

perform_rarefaction_resampling <- function(data_group, n_resamples, sample_size, stat_fn) {
  # Rarefaction logic separated from main function  
}

calculate_confidence_intervals <- function(resampled_values, CI, method = "percentile") {
  # Support different CI methods in future
}
```

## 2. PERFORMANCE IMPROVEMENTS (Medium Priority)

### A. Add Parallel Processing Support
```r
disparity_resample <- function(..., parallel = FALSE, ncores = NULL) {
  # Use safe_parallel_lapply for resampling loops
  # Especially beneficial for large n_resamples
}
```

### B. Memory Optimization Options
```r
disparity_resample <- function(..., return_resampled = TRUE, summary_only = FALSE) {
  # Option to not store all resampled values
  # Just return summary statistics for large datasets
}
```

### C. Progress Indication
```r
disparity_resample <- function(..., progress = TRUE) {
  # Add progress bar for long computations
  # Use txtProgressBar or cli package
}
```

## 3. ENHANCED FUNCTIONALITY (Medium Priority)

### A. Multiple CI Methods
```r
calculate_confidence_intervals <- function(resampled_values, CI, 
                                         method = c("percentile", "bias_corrected", "basic")) {
  # Support bootstrap bias-corrected and accelerated (BCa) intervals
}
```

### B. Custom Statistics Support
```r
disparity_resample <- function(..., custom_statistic = NULL) {
  # Allow users to provide custom statistic functions
  # With validation that they return single numeric values
}
```

### C. Batch Processing
```r
disparity_resample_batch <- function(data_list, ...) {
  # Process multiple datasets with same parameters
  # Return combined results object
}
```

## 4. CODE QUALITY IMPROVEMENTS (High Priority)

### A. Better Error Messages
```r
# Instead of:
stop("CI must be a single numeric value between 0 and 1 (exclusive)")

# Use:
validate_ci <- function(CI) {
  if (!is.numeric(CI) || length(CI) != 1) {
    stop("CI must be a single numeric value, got: ", class(CI)[1], " of length ", length(CI))
  }
  if (CI <= 0 || CI >= 1) {
    stop("CI must be between 0 and 1 (exclusive), got: ", CI)
  }
}
```

### B. Consistent Naming Convention
- Use snake_case consistently throughout
- More descriptive variable names (e.g., `stat_fn` instead of `compute_stat`)

### C. Better Documentation Structure
```r
#' @param Data Input data. See Details for supported formats.
#' @details 
#' ## Supported Data Formats
#' - **Vector**: Univariate analysis (variance only)
#' - **Matrix/Data Frame**: Specimens in rows, variables in columns  
#' - **3D Array**: Landmark data (p × k × n)
#' 
#' ## Available Statistics
#' - `multivariate_variance`: Sum of univariate variances
#' - etc.
```

## 5. TESTING AND ROBUSTNESS (High Priority)

### A. Input Edge Cases
```r
# Add handling for:
# - Single observation groups
# - Identical observations (zero variance)
# - Very small sample sizes
# - Extreme CI values (e.g., 0.99999)
```

### B. Numerical Stability
```r
# Add checks for:
# - Near-singular covariance matrices
# - Overflow/underflow in calculations
# - Degenerate convex hulls
```

## 6. EXAMPLE REFACTORED FUNCTION STRUCTURE

```r
disparity_resample <- function(Data, group = NULL, n_resamples = 1000,
                              statistic = "multivariate_variance", CI = 0.95,
                              bootstrap_rarefaction = "bootstrap", sample_size = NULL,
                              parallel = FALSE, ncores = NULL, progress = TRUE,
                              return_resampled = TRUE, ci_method = "percentile") {
  
  # 1. Input validation (delegated to helper functions)
  inputs <- validate_disparity_inputs(Data, group, n_resamples, statistic, 
                                     CI, bootstrap_rarefaction, sample_size)
  
  # 2. Data preparation (delegated to helper functions)  
  prepared <- prepare_data_and_groups(inputs$Data, inputs$group)
  
  # 3. Get statistic function (delegated to registry)
  stat_fn <- get_statistic_function(statistic, prepared$original_input_is_vector)
  
  # 4. Main resampling loop (simplified)
  results <- perform_resampling_analysis(
    prepared$data, prepared$group_factor, n_resamples, stat_fn,
    bootstrap_rarefaction, sample_size, CI, ci_method, 
    parallel, ncores, progress, return_resampled
  )
  
  # 5. Format output
  format_disparity_results(results, statistic, prepared$single_group_analysis)
}
```

This modular approach would make the code much easier to:
- Test individual components
- Extend with new statistics  
- Maintain and debug
- Read and understand
- Reuse components in other functions
