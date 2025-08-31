#' Resampling-based estimates (bootstrap or rarefaction) of disparity / morphospace occupation
#'
#' Provides a unified interface to obtain resampled (bootstrap or rarefied)
#' estimates of several disparity / morphospace occupation statistics.
#'
#' The function allows choosing among the following multivariate statistics:
#' \itemize{
#'  \item Multivariate variance (trace of the covariance matrix; sum of variances)
#'  \item Mean pairwise Euclidean distance
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
#' Rarefaction requires specifying a valute to the attribute `sample_size`; an error is returned otherwise.
#'
#' @section Observed estimate:
#' For bootstrap resampling the observed estimate reported is the statistic computed on the
#' original (non-resampled) data for each group. For rarefaction, because the purpose is to make
#' groups comparable at a common (smaller) sample size, the observed estimate reported is the
#' mean of the rarefied resampled values for that group (i.e., the mean across all rarefaction
#' replicates).
#'
#' @section Input data types:
#' `Data` can be a data frame, a matrix, a vector, or a 3D array of landmark coordinates
#' (e.g., p landmarks x k dimensions x n specimens). In the latter case, the array is converted
#' internally to a 2D matrix with specimens in rows and (landmark * dimension) variables in columns.
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
#'  confidence intervals.
#' @param bootstrap_rarefaction Either `"bootstrap"` (default) for resampling with replacement or
#'  `"rarefaction"` for resampling without replacement.
#' @param sample_size Either `NULL` (default), a positive integer indicating the number of rows to
#'  sample without replacement when `bootstrap_rarefaction=="rarefaction"`, or the character
#'  `"smallest"` to use the size of the smallest group (all groups rarefied to that size). If
#'  `"smallest"` is supplied but no groups are defined, an error is returned. Required (not `NULL`)
#'  when `bootstrap_rarefaction=="rarefaction"`.
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
#'   # Bootstrap multivariate variance
#'   boot_res = disparity_resample(Data, group=grp, n_resamples=200,
#'                                 statistic="multivariate_variance",
#'                                 bootstrap_rarefaction="bootstrap")
#'   # Direct access to results table
#'   boot_res$results
#'   
#'   # Using the print method for formatted output
#'   print(boot_res)
#'   
#'   # Using the plot method to visualize results
#'   # plot(boot_res)  # Uncomment to create confidence interval plot
#'
#'   # Rarefaction (to the smallest group size) of mean pairwise Euclidean distance
#'   rar_res = disparity_resample(Data, group=grp, n_resamples=200,
#'                                statistic="mean_pairwise_euclidean_distance",
#'                                bootstrap_rarefaction="rarefaction", sample_size="smallest")
#'# Now simulate a third group with larger variance
#'  X3 = MASS::mvrnorm(15, mu=rep(0, 10), Sigma=diag(10)*1.5)
#'  grp2 = factor(c(rep("A", nrow(X1)), rep("B", nrow(X2)), rep("C", nrow(X3))))
#'  boot_res2 = disparity_resample(Data=rbind(X1, X2, X3), group=grp2, n_resamples=1000,
#'                                 statistic="multivariate_variance",
#'                                 bootstrap_rarefaction="bootstrap")
#'  print(boot_res2)
#'  # plot(boot_res2)
#'  # Plot of the obtained (95%) confidence intervals (uncomment to plot)
#' }
#'
#' @import stats
#' @export
disparity_resample=function(Data, group=NULL, n_resamples=1000,
                            statistic="multivariate_variance", CI=0.95,
                            bootstrap_rarefaction="bootstrap", sample_size=NULL) {

  # Input validation ----------------------------------------------------------------
  
  # Validate CI parameter
  if (!is.numeric(CI) || length(CI) != 1 || CI <= 0 || CI >= 1) {
    stop("CI must be a single numeric value between 0 and 1 (exclusive)")
  }
  
  # Validate n_resamples
  if (!is.numeric(n_resamples) || length(n_resamples) != 1 || n_resamples <= 0 || n_resamples != as.integer(n_resamples)) {
    stop("n_resamples must be a positive integer")
  }
  
  # Validate bootstrap_rarefaction
  if (!bootstrap_rarefaction %in% c("bootstrap", "rarefaction")) {
    stop("bootstrap_rarefaction must be either 'bootstrap' or 'rarefaction'")
  }

  # Coerce Data to appropriate matrix / vector form ---------------------------------
  original_input_is_vector=FALSE
  if (is.array(Data) && length(dim(Data))==3) {
    # Assume dimensions: landmarks x dimensions x specimens (p x k x n)
    if (!requireNamespace("Morpho", quietly = TRUE)) {
      stop("Package 'Morpho' is required for 3D array conversion")
    }
    Data=Morpho::vecx(Data)
    # Convert to matrix
  }
  if (is.vector(Data) && !is.list(Data)) {
    Data = as.numeric(Data)
    original_input_is_vector=TRUE
  }
  if (is.data.frame(Data)) { Data = as.matrix(Data) }
  if (!original_input_is_vector && !is.matrix(Data)) {
    stop("Data should be a matrix, data frame, vector, or a 3D array") }

  # Handle missing data -------------------------------------------------------------
  if (original_input_is_vector) {
    missing_indices = which(is.na(Data))
    if (length(missing_indices) > 0) {
      warning(sprintf("Removing %d observations with missing data", length(missing_indices)))
      Data = Data[-missing_indices]
      if (!is.null(group)) {
        group = group[-missing_indices]
      }
    }
  } else {
    missing_rows = which(apply(Data, 1, function(x) any(is.na(x))))
    if (length(missing_rows) > 0) {
      warning(sprintf("Removing %d observations with missing data", length(missing_rows)))
      Data = Data[-missing_rows, , drop=FALSE]
      if (!is.null(group)) {
        group = group[-missing_rows]
      }
    }
  }
  
  # Check for empty data after missing data removal
  n_obs = ifelse(original_input_is_vector, length(Data), nrow(Data))
  if (n_obs == 0) {
    stop("No observations remain after removing missing data")
  }

  # Handle grouping -----------------------------------------------------------------
  n_obs = ifelse(original_input_is_vector, length(Data), nrow(Data))
  
  if (is.null(group)) {
    group_factor = factor(rep("All", n_obs))
    single_group_analysis = TRUE
  } else {
    if (length(group) != n_obs) {
      stop("Length of group does not match number of observations after removing missing data") 
    }
    group_factor = as.factor(group)
    single_group_analysis = FALSE
  }

  group_levels = levels(group_factor)
  
  # Validate group sizes (minimum 2 observations per group for variance calculation)
  group_sizes = table(group_factor)
  min_group_size = min(group_sizes)
  
  if (min_group_size < 2) {
    small_groups = names(group_sizes)[group_sizes < 2]
    stop(sprintf("All groups must have at least 2 observations for variance calculation. Groups with < 2 observations: %s", 
                 paste(small_groups, collapse=", ")))
  }

  # Determine statistic --------------------------------------------------------------
  statistic_choices=c("multivariate_variance", "mean_pairwise_euclidean_distance",
                      "convex_hull_volume", "claramunt_proper_variance")

  if (!original_input_is_vector) {
    if (!(statistic %in% statistic_choices)) {
      stop(paste("statistic should be one of:", paste(statistic_choices, collapse=", "))) 
    }
    
    # Special validation for convex hull volume
    if (statistic == "convex_hull_volume") {
      n_vars = ncol(Data)
      if (min_group_size <= n_vars) {
        stop(sprintf("Convex hull volume requires more observations than variables in each group. Minimum group size: %d, Number of variables: %d", 
                     min_group_size, n_vars))
      }
    }
  }

  stat_label_map=list(multivariate_variance="Multivariate variance",
                      mean_pairwise_euclidean_distance="Mean pairwise Euclidean distance",
                      convex_hull_volume="Convex hull volume",
                      claramunt_proper_variance="Claramunt proper variance",
                      variance_univariate="Variance")

  chosen_statistic = if (original_input_is_vector) {
    stat_label_map$variance_univariate
  } else { stat_label_map[[statistic]] }

  # Helper to compute one statistic -------------------------------------------------
  compute_stat=function(X) {
    if (original_input_is_vector) { return(var(X)) }
    if (statistic=="multivariate_variance") { return(muvar(X)) }
    if (statistic=="mean_pairwise_euclidean_distance") { return(meanpairwiseEuclideanD(X)) }
    if (statistic=="convex_hull_volume") {
      if (!requireNamespace("geometry", quietly = TRUE)) {
        stop("Package 'geometry' is required for convex hull volume") }
      return(geometry::convhulln(X, options="FA")$vol)
    }
  if (statistic=="claramunt_proper_variance") {
      if (!requireNamespace("nlshrink", quietly = TRUE)) {
        stop("Package 'nlshrink' is required for Claramunt proper variance") }
      return(claramunt_proper_variance(X))
    }
  }

  # Rarefaction settings ------------------------------------------------------------
  if (bootstrap_rarefaction=="rarefaction") {
    if (is.null(sample_size)) { 
      stop("sample_size must be provided for rarefaction") 
    }
    if (identical(sample_size, "smallest")) {
      if (length(group_levels)==1) { 
        stop("sample_size='smallest' requires multiple groups") 
      }
      sample_size_num=min(group_sizes)
    } else {
      sample_size_num=as.numeric(sample_size)
      if (is.na(sample_size_num) || sample_size_num<=0 || sample_size_num != as.integer(sample_size_num)) { 
        stop("sample_size must be a positive integer or 'smallest'") 
      }
    }
    
    # Validate that sample_size is appropriate for all groups
    if (sample_size_num < 2) {
      stop("sample_size for rarefaction must be at least 2 for variance calculation")
    }
    
    if (sample_size_num > min_group_size) {
      stop(sprintf("sample_size (%d) is larger than the smallest group size (%d)", 
                   sample_size_num, min_group_size))
    }
    
    # Special validation for convex hull with rarefaction
    if (!original_input_is_vector && statistic == "convex_hull_volume") {
      n_vars = ncol(Data)
      if (sample_size_num <= n_vars) {
        stop(sprintf("Convex hull volume with rarefaction requires sample_size > number of variables. sample_size: %d, Number of variables: %d", 
                     sample_size_num, n_vars))
      }
    }
  }

  # Preallocation of "storage"
  resampled_values_list=list()
  results_rows=list()

  alpha=(1-CI)/2

  for (g in group_levels) {
    idx=which(group_factor==g)
    Xg = if (original_input_is_vector) { Data[idx] } else { Data[idx,,drop=FALSE] }
    n_g = length(idx)

    if (bootstrap_rarefaction=="rarefaction") {
      # sample_size validation already done above
      res_vals=unlist(lapply(seq_len(n_resamples), function(i) {
        sampled_idx=sample(seq_len(n_g), sample_size_num, replace=FALSE)
  if (original_input_is_vector) { compute_stat(Xg[sampled_idx]) } else { compute_stat(Xg[sampled_idx,,drop=FALSE]) }
      }))
      observed_stat=mean(res_vals)
    } else if (bootstrap_rarefaction=="bootstrap") {
      res_vals=unlist(lapply(seq_len(n_resamples), function(i) {
        sampled_idx=sample(seq_len(n_g), n_g, replace=TRUE)
  if (original_input_is_vector) { compute_stat(Xg[sampled_idx]) } else { compute_stat(Xg[sampled_idx,,drop=FALSE]) }
      }))
      observed_stat=compute_stat(Xg)
    } else { stop("bootstrap_rarefaction must be either 'bootstrap' or 'rarefaction'") }

    CI_min=quantile(res_vals, probs=alpha, names=FALSE, type=7)
    CI_max=quantile(res_vals, probs=1-alpha, names=FALSE, type=7)

    resampled_values_list[[g]]=res_vals
    results_rows[[g]]=data.frame(group=g, observed=observed_stat,
                                 CI_min=CI_min, CI_max=CI_max,
                                 row.names=NULL)
  }

  results_df=do.call(rbind, results_rows)
  
  # Fix group handling inconsistency: keep "All" for single group analysis
  if (single_group_analysis) { 
    # Keep "All" as the group label for consistency
    results_df$group = "All"
  }

  # resampled_values output formatting ---------------------------------------------
  resampled_values_out = if (length(group_levels)==1) {
    unname(resampled_values_list[[1]])
  } else { resampled_values_list }

  results_list=list(chosen_statistic=chosen_statistic,
              results=results_df,
              resampled_values=resampled_values_out)
  class(results_list)=c("disparity_resample", "list")
  # Set class for S3 methods
return(results_list)
}


#' Plot method for disparity_resample objects
#'
#' Creates a confidence interval plot for disparity resample results
#'
#' @param x An object of class "disparity_resample"
#' @param ... Additional arguments
#'
#' @return A ggplot object
#' @export
plot.disparity_resample=function(x, ...) {
  # Check if results contain groups or single analysis
  if (any(x$results$group == "All") && nrow(x$results) == 1) {
    # Single group case - already labeled as "All"
    plot_data=x$results
    x_lab="Analysis"
  } else {
    # Multiple groups case
    plot_data=x$results
    x_lab="Group"
  }
  
  # Create plot using internal CI_plot function
  p=CI_plot(data=plot_data, x_var="group", y_var="observed",
            ymin_var="CI_min", ymax_var="CI_max",
            x_lab=x_lab, y_lab=x$chosen_statistic, ...)
  
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
print.disparity_resample=function(x, ...) {
  
  cat("Disparity resampling results\n")
  cat("===========================\n\n")
  cat("Statistic:", x$chosen_statistic, "\n\n")
  
  # Print results table
  print(x$results, row.names=FALSE)
  
  # Check for CI overlap if multiple groups
  if (nrow(x$results) > 1 && !any(x$results$group == "All")) {
    cat("\nConfidence interval overlap assessment:\n")
    
    # Check all pairwise overlaps
    n_groups=nrow(x$results)
    overlaps=matrix(TRUE, nrow=n_groups, ncol=n_groups)
    # Initialize overlap matrix
    
    for (i in seq_len(n_groups)) {
      for (j in seq_len(n_groups)) {
        if (i != j) {
          # Check if CIs overlap: max(min1, min2) <= min(max1, max2)
          overlap_check=max(x$results$CI_min[i], x$results$CI_min[j]) <= 
                        min(x$results$CI_max[i], x$results$CI_max[j])
          overlaps[i, j]=overlap_check
        }
      }
    }
    
    # Summary of overlaps
    all_overlap=all(overlaps[upper.tri(overlaps)])
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


