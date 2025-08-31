#' Resampling-based estimates (bootstrap or rarefaction) of disparity / morphospace occupation
#'
#' Provides a unified interface to obtain resampled (bootstrap or rarefied)
#' estimates of several disparity / morphospace occupation statistics.
#'
#' The function allows choosing among the following multivariate statistics:
#' \itemize{
#'  \item Multivariate variance (trace of the covariance matrix; sum of variances)
#'  \item Mean pairwise Euclidea
#' Plot method for disparity_resample objectstance
#'  \item Convex hull volume (n-dimensional)
#'  \item Claramunt proper variance (Claramunt 2010) based on a linear shrinkage covariance matrix
#' }
#'
#' If the input `Data` is univariate (i.e., a vector), the analysis defaults to computing
#' the (univariate) variance within each group (the attribute `statistic` is ignored).
#'
#' If `bootstrap_rarefaction=="bootstrap"`, the function performs resampling with replacement
#' (i.e., classical bootstrap) by sampling rows of the data matrix / data frame.
#'
#' If `bootstrap_rarefaction=="rarefaction"`, the function performs resampling without replacement
#' at the sample size indicated in `sample_size` (numeric) or, if `sample_size=="smallest"`,
#' at the size of the smallest group (all groups are resampled to that size).
#' Rarefaction requires specifying a value to the attribute `sample_size`; an error is returned otherwise.
#'
#' If `bootstrap_rarefaction=="rarefaction_with_replacement"`, the function performs resampling 
#' with replacement at the sample size indicated in `sample_size`. This is equivalent to 
#' choosing "bootstrap" with a specified `sample_size`.
#'
#' @section Observed estimate:
#' For bootstrap resampling the observed estimate reported is the statistic computed on the
#' original (non-resampled) data for each group. For rarefaction and rarefaction_with_replacement, 
#' because the purpose is to make groups comparable at a common (smaller) sample size, the observed 
#' estimate reported is the mean of the rarefied resampled values for that group (i.e., the mean 
#' across all rarefaction replicates).
#'
#' @section Input data types:
#' `Data` can be a data frame, a matrix, a vector, or a 3D array of landmark coordinates
#' (e.g., p landmarks x k dimensions x n specimens). In the latter case, the array is converted
#' internally to a 2D matrix with specimens in rows and (landmark * dimension) variables in columns.
#'
#' @section Confidence intervals:
#' When `CI = 1`, the function returns the full range (minimum and maximum) of the resampled 
#' values instead of confidence intervals. This can be useful for exploring the complete 
#' variability of the statistic.
#'
#' @note Because of how the computation works, convex hull volume computation requires the number of observations (specimens) to be (substantially) greater than the number of variables (dimensions).
#' In case of shape or similar, consider using the scores along the first (few/several) principal components.
#' Sometimes errors are thrown due to near-zero components, in this case try reducing the number of principal components used.
#' Examples of use of this statistic with geometric morphometric data include Drake & Klingenberg 2010 (American Naturalist), Fruciano et al. 2012 (Environmental Biology of Fishes) and Fruciano et al. 2014 (Biological Journal of the Linnean Society).
#' Because of the sensitivity of this statistic to outliers, usually rarefaction is preferred to bootstrapping.
#'
#' @note "Multivariate variance" is also called "total variance", "Procrustes variance" (in geometric morphometrics) and "sum of univariate variances".
#' Note how the computation here does not divide variance by sample size (other than the normal division performed in the computation of variances).
#'
#' @param Data A data frame, matrix, vector, or 3D array. Observations (specimens) must be in rows
#'  (if a 3D array is supplied, the third dimension is assumed to index specimens).
#' @param group A factor or a vector indicating group membership (will be coerced to factor). If
#'  `NULL` (default) a single analysis is performed on the full dataset.
#' @param n_resamples Number of resampling replicates (default 1000).
#' @param statistic Character string identifying the statistic to compute. One of
#'  `"multivariate_variance"`, `"mean_pairwise_euclidean_distance"`, `"convex_hull_volume"`,
#'  `"claramunt_proper_variance"`. Default is
#'  `"multivariate_variance"`. Ignored for univariate data (vector input).
#' @param CI Desired two-sided confidence interval level (default 0.95) used for percentile
#'  confidence intervals. When `CI = 1`, returns the full range of resampled values.
#' @param bootstrap_rarefaction Either `"bootstrap"` (default) for resampling with replacement, 
#'  `"rarefaction"` for resampling without replacement, or `"rarefaction_with_replacement"` 
#'  for resampling with replacement at a specified sample size.
#' @param sample_size Either `NULL` (default), a positive integer indicating the number of rows to
#'  sample when `bootstrap_rarefaction` is `"rarefaction"` or `"rarefaction_with_replacement"`, 
#'  or the character `"smallest"` to use the size of the smallest group (all groups rarefied to that size). 
#'  If `"smallest"` is supplied but no groups are defined, an error is returned. Required (not `NULL`)
#'  when `bootstrap_rarefaction` is `"rarefaction"` or `"rarefaction_with_replacement"`.
#' @param parallel Logical indicating whether to use parallel processing for resampling (default FALSE).
#' @param ncores Integer number of cores to use for parallel processing. If NULL, defaults to 
#'  `parallel::detectCores() - 1`.
#' @param progress Logical indicating whether to show progress information (default TRUE).
#' @param ci_method Character string specifying confidence interval method. Currently only 
#'  "percentile" is fully supported (default "percentile").
#'
#' @return A list containing:
#'  \describe{
#'    \item{chosen_statistic}{Character vector of length 1 with the human-readable name of the statistic used.}
#'    \item{results}{A data frame with columns `group`, `observed`, `CI_min`, `CI_max`. One row per group.}
#'    \item{resampled_values}{If a single group: numeric vector of length `n_resamples` with the resampled values.
#'      If multiple groups: a named list with one numeric vector (length `n_resamples`) per group.}
#' }
#' 
#' The returned object has class "disparity_resample" and comes with associated S3 methods for 
#' convenient display and visualization:
#' \itemize{
#'   \item \code{\link{print.disparity_resample}}: Prints a formatted summary of results including confidence interval overlap assessment for multiple groups
#'   \item \code{\link{plot.disparity_resample}}: Creates a confidence interval plot using ggplot2
#' }
#'
#' @references Drake AG, Klingenberg CP. 2010. Large-scale diversification of skull shape in domestic dogs: disparity and modularity. American Naturalist 175:289-301.
#' @references Claramunt S. 2010. Discovering exceptional diversifications at continental scales: the case of the endemic families of Neotropical Suboscine passerines. Evolution 64:2004-2019.
#' @references Fruciano C, Tigano C, Ferrito V. 2012. Body shape variation and colour change during growth in a protogynous fish. Environmental Biology of Fishes 94:615-622.
#' @references Fruciano C, Pappalardo AM, Tigano C, Ferrito V. 2014. Phylogeographical relationships of Sicilian brown trout and the effects of genetic introgression on morphospace occupation. Biological Journal of the Linnean Society 112:387-398.
#' @references Fruciano C, Franchini P, Raffini F, Fan S, Meyer A. 2016. Are sympatrically speciating Midas cichlid fish special? Patterns of morphological and genetic variation in the closely related species Archocentrus centrarchus. Ecology and Evolution 6:4102-4114.
#' @seealso \code{\link{disparity_test}}, \code{\link{print.disparity_resample}}, \code{\link{plot.disparity_resample}}
#'
#' @examples
#' set.seed(123)
#' # Simulate two groups with different means but same covariance
#' if (requireNamespace("MASS", quietly = TRUE)) {
#'   X1 = MASS::mvrnorm(20, mu=rep(0, 10), Sigma=diag(10))
#'   X2 = MASS::mvrnorm(35, mu=rep(2, 10), Sigma=diag(10))
#'   Data = rbind(X1, X2)
#'   grp = factor(c(rep("A", nrow(X1)), rep("B", nrow(X2))))
#'
#'   # Bootstrap multivariate variance with BCa intervals (if boot package available)
#'   boot_res = disparity_resample(Data, group=grp, n_resamples=200,
#'                                 statistic="multivariate_variance",
#'                                 bootstrap_rarefaction="bootstrap",
#'                                 ci_method="bca")
#'   print(boot_res)
#'   
#'   # Rarefaction (to the smallest group size) of mean pairwise Euclidean distance
#'   rar_res = disparity_resample(Data, group=grp, n_resamples=200,
#'                                statistic="mean_pairwise_euclidean_distance",
#'                                bootstrap_rarefaction="rarefaction", 
#'                                sample_size="smallest")
#'   
#'   # Rarefaction with replacement at fixed sample size
#'   raref_repl_res = disparity_resample(Data, group=grp, n_resamples=200,
#'                                       statistic="multivariate_variance",
#'                                       bootstrap_rarefaction="rarefaction_with_replacement",
#'                                       sample_size=15)
#'   
#'   # Get full range of values (CI = 1)
#'   range_res = disparity_resample(Data, group=grp, n_resamples=200,
#'                                  statistic="multivariate_variance",
#'                                  bootstrap_rarefaction="bootstrap",
#'                                  CI=1.0)
#'   
#'   # Parallel processing with multiple cores
#'   # parallel_res = disparity_resample(Data, group=grp, n_resamples=1000,
#'   #                                   statistic="multivariate_variance",
#'   #                                   bootstrap_rarefaction="bootstrap",
#'   #                                   ncores=4)
#' }
#'
#' @import stats
#' @export
disparity_resample = function(Data, group = NULL, n_resamples = 1000,
                             statistic = "multivariate_variance", CI = 0.95,
                             bootstrap_rarefaction = "bootstrap", sample_size = NULL,
                             parallel = FALSE, ncores = NULL, progress = TRUE,
                             ci_method = "percentile") {

  # Input validation ----------------------------------------------------------------
  inputs = validate_disparity_inputs(Data, group, n_resamples, statistic, CI, 
                                    bootstrap_rarefaction, sample_size, parallel, ncores)
  
  # Data preparation ----------------------------------------------------------------
  prepared = prepare_data_and_groups(inputs$Data, inputs$group)
  
  # Set up parallel processing
  if (inputs$parallel && is.null(inputs$ncores)) {
    inputs$ncores = max(1, parallel::detectCores(logical = FALSE) - 1)
  }
  
  # Group size validation
  validate_group_sizes(prepared$group_factor, inputs$statistic, 
                      prepared$original_input_is_vector, 
                      if (!prepared$original_input_is_vector) prepared$data else NULL)

  # Determine statistic function
  stat_fn = get_statistic_function(inputs$statistic, prepared$original_input_is_vector)
  registry = create_statistic_registry()
  chosen_statistic = if (prepared$original_input_is_vector) {
    registry$variance_univariate$name
  } else { 
    registry[[inputs$statistic]]$name 
  }

  # Handle rarefaction settings -----------------------------------------------------
  sample_size_num = NULL
  if (inputs$bootstrap_rarefaction %in% c("rarefaction", "rarefaction_with_replacement")) {
    if (is.null(inputs$sample_size)) { 
      stop("sample_size must be provided when bootstrap_rarefaction is '", 
           inputs$bootstrap_rarefaction, "'") 
    }
    
    group_levels = levels(prepared$group_factor)
    group_sizes = table(prepared$group_factor)
    
    if (identical(inputs$sample_size, "smallest")) {
      if (length(group_levels) == 1) { 
        stop("sample_size='smallest' requires multiple groups") 
      }
      sample_size_num = min(group_sizes)
    } else {
      sample_size_num = as.numeric(inputs$sample_size)
      if (is.na(sample_size_num) || sample_size_num <= 0 || 
          sample_size_num != as.integer(sample_size_num)) { 
        stop("sample_size must be a positive integer or 'smallest'") 
      }
      sample_size_num = as.integer(sample_size_num)
    }
    
    # Re-validate group sizes with sample_size
    validate_group_sizes(prepared$group_factor, inputs$statistic, 
                        prepared$original_input_is_vector, 
                        if (!prepared$original_input_is_vector) prepared$data else NULL,
                        sample_size_num)
  }

  # Main resampling loop ------------------------------------------------------------
  group_levels = levels(prepared$group_factor)
  resampled_values_list = list()
  results_rows = list()

  if (inputs$progress && length(group_levels) > 1) {
    cat("Processing", length(group_levels), "groups...\n")
  }

  for (i in seq_along(group_levels)) {
    g = group_levels[i]
    if (inputs$progress && length(group_levels) > 1) {
      cat("Group", i, "of", length(group_levels), ":", g, "\n")
    }
    
    idx = which(prepared$group_factor == g)
    Xg = if (prepared$original_input_is_vector) { 
      prepared$data[idx] 
    } else { 
      prepared$data[idx, , drop = FALSE] 
    }

    # Perform resampling based on method
    if (inputs$bootstrap_rarefaction == "bootstrap") {
      res_vals = perform_bootstrap_resampling(
        Xg, inputs$n_resamples, stat_fn, prepared$original_input_is_vector,
        sample_size = NULL, parallel = inputs$parallel, ncores = inputs$ncores
      )
      observed_stat = stat_fn(Xg)
      
    } else if (inputs$bootstrap_rarefaction == "rarefaction") {
      res_vals = perform_rarefaction_resampling(
        Xg, inputs$n_resamples, sample_size_num, stat_fn, prepared$original_input_is_vector,
        parallel = inputs$parallel, ncores = inputs$ncores
      )
      observed_stat = mean(res_vals)
      
    } else if (inputs$bootstrap_rarefaction == "rarefaction_with_replacement") {
      res_vals = perform_bootstrap_resampling(
        Xg, inputs$n_resamples, stat_fn, prepared$original_input_is_vector,
        sample_size = sample_size_num, parallel = inputs$parallel, ncores = inputs$ncores
      )
      observed_stat = mean(res_vals)
      
    } else { 
      stop("Unknown bootstrap_rarefaction method: ", inputs$bootstrap_rarefaction) 
    }

    # Calculate confidence intervals
    ci_vals = calculate_confidence_intervals(res_vals, inputs$CI, ci_method)

    resampled_values_list[[g]] = res_vals
    results_rows[[g]] = data.frame(
      group = g, 
      observed = observed_stat,
      CI_min = ci_vals["CI_min"], 
      CI_max = ci_vals["CI_max"],
      row.names = NULL
    )
  }

  # Format results ------------------------------------------------------------------
  results_df = do.call(rbind, results_rows)
  
  # Keep "All" for single group analysis for consistency
  if (prepared$single_group_analysis) { 
    results_df$group = "All"
    names(resampled_values_list) = "All"
  }

  # Format resampled values output
  resampled_values_out = if (length(group_levels) == 1) {
    unname(resampled_values_list[[1]])
  } else { 
    resampled_values_list 
  }

  # Create results object
  results_list = list(
    chosen_statistic = chosen_statistic,
    results = results_df,
    resampled_values = resampled_values_out
  )
  class(results_list) = c("disparity_resample", "list")
  
  if (inputs$progress) {
    cat("Analysis complete.\n")
  }
  
  return(results_list)
}


#' Plot method for disparity_resample objects
#'
#' Creates a confidence interval plot for disparity resample results
#'
#' @param x An object of class "disparity_resample"
#' @param point_color A single color or vector of colors for point estimates passed to internal CI_plot.
#'   If a vector, rules follow CI_plot: length 1 (recycled), length equal to number of x levels, or number of rows.
#'   (default "darkblue").
#' @param errorbar_color A single color or vector of colors for error bars (same length rules as point_color; default "darkred").
#' @param ... Additional arguments passed to CI_plot
#'
#' @return A ggplot object
#' @export
plot.disparity_resample = function(x, point_color = "darkblue", errorbar_color = "darkred", ...) {
  # Check if results contain groups or single analysis
  if (any(x$results$group == "All") && nrow(x$results) == 1) {
    # Single group case - already labeled as "All"
    plot_data = x$results
    x_lab = "Analysis"
  } else {
    # Multiple groups case
    plot_data = x$results
    x_lab = "Group"
  }
  
  # Create plot using internal CI_plot function
  p = CI_plot(data = plot_data, x_var = "group", y_var = "observed",
            ymin_var = "CI_min", ymax_var = "CI_max",
            x_lab = x_lab, y_lab = x$chosen_statistic,
            point_color = point_color, errorbar_color = errorbar_color, ...)
  
  return(p)
}


#' Print method for disparity_resample objects
#'
#' Prints results table and checks for CI overlap among groups
#'
#' @param x An object of class "disparity_resample"
#' @param ... Additional arguments (not used)
#'
#' @return Invisibly returns the input object
#' @export
print.disparity_resample = function(x, ...) {
  
  cat("Disparity resampling results\n")
  cat("===========================\n\n")
  cat("Statistic:", x$chosen_statistic, "\n\n")
  
  # Print results table
  print(x$results, row.names = FALSE)
  
  # Check for CI overlap if multiple groups
  if (nrow(x$results) > 1 && !any(x$results$group == "All")) {
    cat("\nConfidence interval overlap assessment:\n")
    
    # Check all pairwise overlaps
    n_groups = nrow(x$results)
    overlaps = matrix(TRUE, nrow = n_groups, ncol = n_groups)
    # Initialize overlap matrix
    
    for (i in seq_len(n_groups)) {
      for (j in seq_len(n_groups)) {
        if (i != j) {
          # Check if CIs overlap: max(min1, min2) <= min(max1, max2)
          overlap_check = max(x$results$CI_min[i], x$results$CI_min[j]) <= 
                        min(x$results$CI_max[i], x$results$CI_max[j])
          overlaps[i, j] = overlap_check
        }
      }
    }
    
    # Summary of overlaps
    all_overlap = all(overlaps[upper.tri(overlaps)])
    # Check upper triangle only (avoid diagonal and duplicates)
    
    if (all_overlap) {
      cat("All confidence intervals overlap.\n")
    } else {
      cat("At least one pair of confidence intervals does not overlap.\n")
    }
  }
  
  cat("\n")
  invisible(x)
}