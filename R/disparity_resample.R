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
#' @references Drake AG, Klingenberg CP. 2010. Large-scale diversification of skull shape in domestic dogs: disparity and modularity. American Naturalist 175:289-301.
#' @references Claramunt S. 2010. Discovering exceptional diversifications at continental scales: the case of the endemic families of Neotropical Suboscine passerines. Evolution 64:2004-2019.
#' @references Fruciano C, Pappalardo AM, Tigano C, Ferrito V. 2014. Phylogeographical relationships of Sicilian brown trout and the effects of genetic introgression on morphospace occupation. Biological Journal of the Linnean Society 112:387-398.
#'
#' @seealso \code{\link{disparity_test}}
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
#'   boot_res$results
#'
#'   # Rarefaction (to the smallest group size) of mean pairwise Euclidean distance
#'   rar_res = disparity_resample(Data, group=grp, n_resamples=200,
#'                                statistic="mean_pairwise_euclidean_distance",
#'                                bootstrap_rarefaction="rarefaction", sample_size="smallest")
#'   rar_res$results
#' }
#'
#' @import stats
#' @export
disparity_resample=function(Data, group=NULL, n_resamples=1000,
                            statistic="multivariate_variance", CI=0.95,
                            bootstrap_rarefaction="bootstrap", sample_size=NULL) {

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

  # Handle grouping -----------------------------------------------------------------
  if (is.null(group)) {
    group_factor = factor(rep("All", ifelse(original_input_is_vector, length(Data), nrow(Data))))
  } else {
    if (length(group)!=(ifelse(original_input_is_vector, length(Data), nrow(Data)))) {
      stop("Length of group does not match number of observations") }
    group_factor = as.factor(group)
  }

  group_levels = levels(group_factor)

  # Determine statistic --------------------------------------------------------------
  statistic_choices=c("multivariate_variance", "mean_pairwise_euclidean_distance",
                      "convex_hull_volume", "claramunt_proper_variance")

  if (!original_input_is_vector) {
    if (!(statistic %in% statistic_choices)) {
      stop(paste("statistic should be one of:", paste(statistic_choices, collapse=", "))) }
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
    if (is.null(sample_size)) { stop("sample_size must be provided for rarefaction") }
    if (identical(sample_size, "smallest")) {
      if (length(group_levels)==1) { stop("sample_size='smallest' requires multiple groups") }
      group_sizes=table(group_factor)
      sample_size_num=min(group_sizes)
    } else {
      sample_size_num=as.numeric(sample_size)
  if (is.na(sample_size_num) || sample_size_num<=0) { stop("sample_size must be a positive integer or 'smallest'") }
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
      if (sample_size_num>n_g) { stop("sample_size is larger than group sample size for group ", g) }
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
  if (is.null(group)) { results_df$group=NA }

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
  if (all(is.na(x$results$group))) {
    # Single group case - create simple label
    plot_data=x$results
    plot_data$group="All"
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
  if (nrow(x$results) > 1 && !all(is.na(x$results$group))) {
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


