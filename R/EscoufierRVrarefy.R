#' Rarefied version of Escoufier RV coefficient
#'
#' Computes a rarefied estimate of the Escoufier RV coefficient
#' to account for the dependence on sample size
#' of RV, so that RV can be compared among groups with
#' different number of observations (sample sizes)
#'
#' This function computes a rarefied estimate of Escoufier RV coefficient
#' as suggested by Fruciano et al 2013 - Plos One
#' This can be useful to compare RV among groups with the same variables
#' but different sample sizes (as RV depends on sample size, see Fruciano et al 2013,
#' where this procedure is described).
#' The idea is the one rarefies the two groups at the same sample size.
#' Following the approach in Fruciano et al. 2013 - Plos One, "rarefaction" 
#' is meant resampling to a smaller sample size with replacement (like in bootstrap).
#'
#' @section Notice:
#' the function does NOT perform GPA on each rarefied sample
#' this may or may not make a difference in estimates.
#' In many cases, it will probably not make much difference
#' (e.g., Fig. 2 in Fruciano et al 2013)
#'
#' @section Citation:
#' If you use this function please cite both
#' Fruciano et al. 2013 (for using the rarefaction procedure)
#' and Escoufier 1973 (because the procedure is based on Escoufier RV)
#'
#' @param Block1,Block2 Matrices or data frames
#' containing each block of variables
#' (observations in rows, variables in columns).
#' @param reps number of resamplings to obtain the rarefied estimate
#' @param samplesize sample size to which the rarefaction procedure is carried out
#' @param group factor or vector indicating group membership for observations.
#' If NULL (default), a single-group analysis is performed. If provided,
#' the analysis is performed separately for each group.
#' @param CI confidence interval level (default 0.95). Must be between 0 and 1.
#'
#' @seealso \code{\link{EscoufierRV}}
#' @return The function outputs a list with the following elements:
#'  \describe{
#'   \item{results}{A data frame with columns `group`, `Mean`, `Median`, `CI_min`, `CI_max`. 
#'                  One row per group. For single-group analysis, group will be "All".}
#'   \item{AllRarefiedSamples}{A list with all RV values obtained using the rarefaction procedure for each group.
#'                             For single-group analysis, this will be a list with one element named "All".}
#' }
#' 
#' The returned object has class "EscoufierRVrarefy" and comes with associated S3 methods for 
#' convenient display and visualization:
#' \itemize{
#'   \item \code{\link{print.EscoufierRVrarefy}}: Prints a formatted summary of results including confidence interval overlap assessment for multiple groups
#'   \item \code{\link{plot.EscoufierRVrarefy}}: Creates a confidence interval plot using ggplot2
#' }
#'
#' @references Escoufier Y. 1973. Le Traitement des Variables Vectorielles. Biometrics 29:751-760.
#' @references Fruciano C, Franchini P, Meyer A. 2013. Resampling-Based Approaches to Study Variation in Morphological Modularity. PLoS ONE 8:e69376.
#' @examples
#' library(MASS)
#' set.seed(123)
#' Pop=mvrnorm(100000,mu=rep(0,100), Sigma=diag(100))
#' # Create a population of 100,000 'individuals'
#' # as multivariate normal random data
#' # We will consider the first 20 columns as the first
#' # block of variables, and the following one as the second block
#'
#' A=Pop[1:50,]
#' B=Pop[501:700,]
#' # Take two groups (A and B)
#' # from the same population (there should be no difference
#' # between them)
#'
#' EscoufierRV(A[,1:20],A[,21:ncol(A)])
#' EscoufierRV(B[,1:20],B[,21:ncol(B)])
#' # Notice how we obtain very different values of Escoufier RV
#' # (this is because they two groups have very different
#' # sample sizes, one 50 observations, the other 200)
#'
#' RarA=RVrarefied(A[,1:20],A[,21:ncol(A)],reps=1000,samplesize=30)
#' RarB=RVrarefied(B[,1:20],B[,21:ncol(A)],reps=1000,samplesize=30)
#' RarA$results  # Data frame with Mean, Median, CI_min, CI_max
#' RarB$results  # Data frame with Mean, Median, CI_min, CI_max
#' # Rarefying both groups at the same sample size
#' # (in this case 30)
#' # it is clear that the two groups have very similar levels
#' # of association between blocks
#'
#' # Multi-group analysis with custom CI
#' combined_data = rbind(A, B)
#' group_labels = c(rep("GroupA", nrow(A)), rep("GroupB", nrow(B)))
#' multi_result = RVrarefied(combined_data[,1:20], combined_data[,21:ncol(combined_data)], 
#'                          reps=1000, samplesize=30, group=group_labels, CI=0.90)
#' print(multi_result$results)  # Data frame with results for each group
#' # Columns: group, Mean, Median, CI_min, CI_max
#'
#' @export
RVrarefied = function(Block1, Block2, reps = 1000, samplesize, group = NULL, CI = 0.95) {
  # Basic input validation and coercion
  validated = validate_blocks(Block1, Block2, samplesize)
  Block1 = validated$Block1
  Block2 = validated$Block2

  if (!is.null(group) && length(group) != nrow(Block1)) {
    stop("Error: group vector must have the same length as the number of observations")
  }

  if (CI <= 0 || CI >= 1) {
    stop("Error: CI must be between 0 and 1")
  }

  alpha = (1 - CI) / 2  # For two-sided CI

  # Handle single group analysis by treating as multi-group with one group
  if (is.null(group)) {
    group = rep("All", nrow(Block1))  # Create single group named "All"
  }
    
    # Multi-group analysis (unified approach)
    group = as.factor(group)
    group_levels = levels(group)
    results_rows = list()
    resampled_values_list = list()
    
    for (g in group_levels) {
        group_indices = which(group == g)
        Block1_g = Block1[group_indices, , drop = FALSE]
        Block2_g = Block2[group_indices, , drop = FALSE]
        
        if (nrow(Block1_g) < samplesize) {
            stop(paste("Error: group", g, "has fewer observations than the requested sample size"))
        }
    # Pre-allocate and sample by row indices to avoid repeated cbind
    n_g = nrow(Block1_g)
    RV_g = vector(length = reps)

    for (i in 1:reps) {
      idx = sample.int(n_g, samplesize, replace = TRUE)
      RV_g[i] = EscoufierRV(Block1_g[idx, , drop = FALSE], Block2_g[idx, , drop = FALSE])
    }
        
        if(TRUE %in% is.na(RV_g)) {
        warning(paste("Some of the rarefied RV for group", g, "have given NA \n
            This might be due to random samples all with the same observations\n
            Please, inspect your data and how frequent it is in the RV computed for rarefied samples\n
            to ensure that the results are still meaningful
            "))
        }
        
        RV_g_clean = na.omit(RV_g)
        # remove any NAs
        
        # Calculate statistics
        mean_rv = mean(RV_g_clean)
        median_rv = median(RV_g_clean)
        CI_min = quantile(RV_g_clean, probs = alpha, names = FALSE, type = 7)
        CI_max = quantile(RV_g_clean, probs = 1 - alpha, names = FALSE, type = 7)
        
        # Store results
        results_rows[[g]] = data.frame(group = g, Mean = mean_rv, Median = median_rv,
                                     CI_min = CI_min, CI_max = CI_max,
                                     row.names = NULL)
        resampled_values_list[[g]] = RV_g
    }
    
    results_df = do.call(rbind, results_rows)
    
  Results = list(results = results_df, AllRarefiedSamples = resampled_values_list)
  class(Results) = c("EscoufierRVrarefy", "list")
  attr(Results, "CI") = CI
  return(Results)
}


# Helper: validate and coerce input blocks
validate_blocks = function(Block1, Block2, samplesize) {
  # Ensure matrix-like and numeric
  if (!is.matrix(Block1)) Block1 = as.matrix(Block1)
  if (!is.matrix(Block2)) Block2 = as.matrix(Block2)

  if (!is.numeric(Block1)) storage.mode(Block1) = "double"
  if (!is.numeric(Block2)) storage.mode(Block2) = "double"

  if (nrow(Block1) != nrow(Block2)) {
    stop("Error: the two blocks should have the same number of rows (observations)")
  }

  if (!is.numeric(samplesize) || length(samplesize) != 1 || samplesize < 1) {
    stop("Error: samplesize must be a single numeric value >= 1")
  }

  # Check for finite values
  if (any(!is.finite(Block1))) stop("Error: Block1 contains non-finite values (NA/NaN/Inf)")
  if (any(!is.finite(Block2))) stop("Error: Block2 contains non-finite values (NA/NaN/Inf)")

  # Check for constant columns (zero variance)
  var1 = apply(Block1, 2, stats::var)
  var2 = apply(Block2, 2, stats::var)
  if (any(var1 == 0)) warning("Warning: Block1 contains column(s) with zero variance")
  if (any(var2 == 0)) warning("Warning: Block2 contains column(s) with zero variance")

  return(list(Block1 = Block1, Block2 = Block2))
}



#' Escoufier RV coefficient
#'
#' Computes the Escoufier RV coefficient
#'
#'
#' This function computes the usual version of the Escoufier RV coefficient (Escoufier, 1973),
#' which quantifies the level of association between two multivariate blocks
#' of variables. The function accepts two blocks of variables, either two data frames
#' or two matrices each of n observations (specimens) as rows.
#' The two blocks must have the same number of rows (specimens), but can have
#' different number of columns (variables, such as landmark coordinates).
#' The Escoufier RV has been shown (Fruciano et al. 2013) to be affected
#' by sample size so comparisons of groups (e.g., species, populations)
#' with different sample size should be avoided, unless steps are taken to account
#' for this problem
#'
#' @param Block1,Block2 Matrices or data frames
#' containing each block of variables
#' (observations in rows, variables in columns).
#'
#' @return The function returns a number, corresponding to the Escoufier RV coefficient
#' @references Escoufier Y. 1973. Le Traitement des Variables Vectorielles. Biometrics 29:751-760.
#' @references Fruciano C, Franchini P, Meyer A. 2013. Resampling-Based Approaches to Study Variation in Morphological Modularity. PLoS ONE 8:e69376.
#'
#' @seealso \code{\link{RVrarefied}}
#'
#' @examples
#' library(MASS)
#' set.seed(123)
#' A=mvrnorm(100,mu=rep(0,100), Sigma=diag(100))
#' # Create a sample of 100 'individuals'
#' # as multivariate normal random data
#' # We will consider the first 20 columns as the first
#' # block of variables, and the following one as the second block
#'
#' EscoufierRV(A[,1:20],A[,21:ncol(A)])
#' # Compute the EscoufierRV using the two blocks of variables
#'
#' @export
EscoufierRV = function(Block1, Block2) {
    if (nrow(Block1) != nrow(Block2)) {
        stop(paste("Error: the two blocks should have the same number of rows (observations)"))
    } else {
      tS1S1=sum(cov(Block1)^2)
      tS2S2=sum(cov(Block2)^2)
      S1S2=cov(Block1, Block2)
      RV=sum(S1S2^2)/sqrt(tS1S1*tS2S2)
        return(RV)
    }
}


####################
#### S3 methods ####
####################


#' Plot method for EscoufierRVrarefy objects
#'
#' Creates a confidence interval plot for RV rarefaction results
#'
#' @param x An object of class "EscoufierRVrarefy"
#' @param point_color A single color or a vector of colors for point estimates.
#'   If length 1, the same color is used for all points. If length equals the
#'   number of groups, colors are assigned per group. (default "darkblue")
#' @param errorbar_color A single color or a vector of colors for error bars.
#'   Follows the same recycling rules as `point_color`. (default "darkred")
#' @param ... Additional arguments passed to the underlying plotting function
#'
#' @return A ggplot object
#' @export
plot.EscoufierRVrarefy = function(x, point_color = "darkblue", errorbar_color = "darkred", ...) {
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
  p = CI_plot(data = plot_data, x_var = "group", y_var = "Mean",
              ymin_var = "CI_min", ymax_var = "CI_max",
              x_lab = x_lab, y_lab = "Rarefied RV", 
              point_color = point_color, errorbar_color = errorbar_color, ...)
  
  return(p)
}


#' Print method for EscoufierRVrarefy objects
#'
#' Prints results table and checks for CI overlap among groups
#'
#' @param x An object of class "EscoufierRVrarefy"
#' @param ... Additional arguments (not used)
#'
#' @return Invisibly returns the input object
#' @export
print.EscoufierRVrarefy = function(x, ...) {
  
  cat("Rarefied Escoufier RV coefficient results\n")
  cat("=========================================\n\n")
  cat("Statistic: Rarefied Escoufier RV coefficient\n")
  
  # Note about CI interpretation - read from attribute if present
  CI_attr = attr(x, "CI")
  if (!is.null(CI_attr)) {
    cat(paste0("Confidence level: ", format(100 * CI_attr, digits = 4), "%\n\n"))
  } else {
    cat("Confidence level: (not stored)\n\n")
  }
  
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
  
  # Report NA counts per group if present in AllRarefiedSamples
  if (!is.null(x$AllRarefiedSamples)) {
    na_counts = vapply(x$AllRarefiedSamples, function(v) sum(is.na(v)),
                       FUN.VALUE = integer(1))
    if (any(na_counts > 0)) {
      cat('\nNumber of NA rarefied RV values per group:\n')
      print(data.frame(group = names(na_counts),
                       NA_count = as.integer(na_counts)),
            row.names = FALSE)
    }
  }
  
  cat("\n")
  invisible(x)
}
