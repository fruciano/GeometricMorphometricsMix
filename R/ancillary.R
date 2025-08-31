### Ancillary functions
### (not exported)


# Functions to repeat rows or columns
#' @noRd
rep_row=function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
#' @noRd
rep_col=function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}


# Functions to convert radians to degrees
# and degrees to radians
rad2deg = function(rad) {(rad * 180) / (pi)}
deg2rad = function(deg) {(deg * pi) / (180)}


# Function to compute "multivariate variance"
# (sum of univariate variances)
muvar=function(X) {
  sum(apply(X, 2, var))
}

# Function to compute "Proper variance" as described in
# Claramunt 2010 - Evolution
# It uses a linear shrinkage estimator of the covariance matrix
claramunt_proper_variance=function(X) {
  covX=nlshrink::linshrink_cov(as.matrix(X))
  eigenv=eigen(covX)$values
  proper_variance=(sum(sqrt(eigenv)))^2
return(proper_variance)
}


# Function to compute mean pairwise Euclidean distances
meanpairwiseEuclideanD=function(X) {
  mean(dist(X,method = "euclidean"))
}

# Function to compute the Euclidean distance
# between two observations
pdistance=function(X1, X2) {
  dist(rbind(X1, X2),method = "euclidean")
}


#' Cross-platform safe parallel lapply
#'
#' Select an appropriate parallel backend depending on the platform and run a
#' parallelized lapply safely. On Windows a PSOCK cluster is used
#' (parallel::parLapply). On Unix-like systems the function will try
#' to use forking via parallel::mclapply and will fall back to a PSOCK
#' cluster if forking is not available or fails. If \code{ncores == 1}
#' this simply calls \code{lapply}. 
#'
#' This function is used internally by a few functions in this package but is exported
#' so that users can reuse it for lightweight cross-platform parallel
#' work (for example on Windows where forking is not available).
#'
#' @param X A list or vector to iterate over.
#' @param FUN A function to apply to each element of \code{X}.
#' @param ncores Integer number of worker processes to use. Defaults to
#'   \code{(parallel::detectCores(logical = FALSE)-1)} and is limited to
#'   \code{length(X)}.
#' @param packages Character vector of package names to be required on the
#'   worker processes (each worker will try to \code{require()} them).
#' @param export Character vector of names of objects in the current
#'   environment to export to worker environments via \code{clusterExport}.
#' @param seed Optional integer seed for RNG on workers. When using PSOCK
#'   clusters this is passed to \code{clusterSetRNGStream}; with forking
#'   (mclapply) the seed is set via \code{set.seed()} on the master prior to
#'   forking.
#' @param type One of \code{"auto"}, \code{"mclapply"} or \code{"psock"}
#'   to force the backend selection. Default is \code{"auto"}.
#' @return A \code{list} with the results of applying \code{FUN} to each
#'   element of \code{X} (same shape as \code{lapply}).
#' @examples
#' \dontrun{
#' # simple usage
#' safe_parallel_lapply(1:4, function(i) i^2, ncores = 2)
#'
#' # more complex example: simulate work with a pause to compare serial vs parallel timing
#' slow_task = function(i) {
#'   # simulate a time-consuming operation (0.5 seconds)
#'   Sys.sleep(0.5)
#'   return(i^2)
#' }
#'
#' # run serially (ncores = 1)
#' t_serial = system.time({
#'   res_serial = safe_parallel_lapply(1:8, slow_task, ncores = 1)
#' })
#'
#' # run in parallel (ncores = 4)
#' t_parallel = system.time({
#'   res_parallel = safe_parallel_lapply(1:8, slow_task, ncores = 4)
#' })
#'
#' # show timings
#' print(t_serial)    # expected ~ 4 seconds (8 * 0.5)
#' print(t_parallel)  # expected substantially less on a 4-core machine
#'
#' # verify results are identical
#' identical(unlist(res_serial), unlist(res_parallel))
#' }
#' @export
safe_parallel_lapply = function(X, FUN, ncores = (parallel::detectCores(logical = FALSE)-1),
                                packages = NULL, export = NULL, seed = NULL,
                                type = c("auto", "mclapply", "psock")) {
  type = match.arg(type)
  ncores = as.integer(max(1, min(length(X), ncores)))

  # basic input validation
  if (!is.function(FUN)) stop("'FUN' must be a function")
  if (!is.null(export)) export = as.character(export)
  if (!is.null(packages)) packages = as.character(packages)

  if (ncores <= 1) return(lapply(X, FUN))

  os_type = .Platform$OS.type  # "windows" or "unix"

  # Helper to run on PSOCK cluster
  run_psock = function() {
    cl = parallel::makeCluster(ncores, type = "PSOCK")
    # ensure cluster is stopped on exit of this function
    on.exit(parallel::stopCluster(cl), add = TRUE)

    if (!is.null(export)) parallel::clusterExport(cl, export, envir = parent.frame())
    if (!is.null(packages)) {
      # export the packages vector so workers can access it
      parallel::clusterExport(cl, varlist = "packages", envir = parent.frame())
      parallel::clusterEvalQ(cl, {
        for (pk in packages) {
          try(require(pk, character.only = TRUE), silent = TRUE)
        }
        NULL
      })
    }
    if (!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)

    # run in a tryCatch to ensure cluster is stopped and errors are informative
    res = tryCatch({
      parallel::parLapply(cl, X, FUN)
    }, error = function(e) {
      stop("Error while running parLapply on PSOCK cluster: ", conditionMessage(e))
    })
    return(res)
  }

  # If user forced a backend
  if (type == "mclapply") {
    if (os_type == "windows") {
      # cannot fork on windows; fallback to psock
      return(run_psock())
    }
    if (!is.null(seed)) set.seed(seed)
    res = parallel::mclapply(X, FUN, mc.cores = ncores, mc.preschedule = TRUE, mc.set.seed = TRUE)
    return(res)
  }

  if (type == "psock") return(run_psock())

  # auto: choose best available
  if (os_type == "windows") {
    return(run_psock())
  } else {
    # try mclapply but fallback to PSOCK on error
    try_res = try({
      if (!is.null(seed)) set.seed(seed)
      parallel::mclapply(X, FUN, mc.cores = ncores, mc.preschedule = TRUE, mc.set.seed = TRUE)
    }, silent = TRUE)
    if (!inherits(try_res, "try-error")) return(try_res)
    # fallback
    return(run_psock())
  }
}


# General resampling functions ===========================================

#' Perform unified resampling with confidence intervals
#'
#' Unified wrapper function that performs resampling using different methods
#' and calculates confidence intervals. For bootstrap, uses internal resampling
#' functions to avoid Windows parallel issues, but can create boot objects
#' for bias-corrected accelerated (BCa) intervals when requested.
#' For rarefaction methods, uses internal resampling functions.
#'
#' @param data_group Vector or matrix of data for resampling
#' @param n_resamples Number of resampling replicates
#' @param stat_fn Function to compute on each resampled dataset
#' @param original_input_is_vector Logical indicating if input was a vector
#' @param bootstrap_rarefaction Character: "bootstrap", "rarefaction", or "rarefaction_with_replacement"
#' @param sample_size Optional integer for fixed sample size (required for rarefaction methods)
#' @param CI Confidence interval level (0 < CI <= 1)
#' @param ci_method Character string specifying CI method ("percentile", "bca" for bootstrap only)
#' @param parallel Logical indicating whether to use parallel processing (default: FALSE)
#' @param ncores Number of cores for parallel processing (ignored if parallel=FALSE)
#' @return List with elements: resampled_values, observed_stat, CI_min, CI_max
#' @examples
#' \dontrun{
#'   # Bootstrap with BCa confidence intervals
#'   x <- rnorm(50)
#'   result <- perform_unified_resampling(x, 1000, mean, TRUE, "bootstrap", 
#'                                       CI = 0.95, ci_method = "bca")
#'   
#'   # Rarefaction
#'   result2 <- perform_unified_resampling(x, 1000, mean, TRUE, "rarefaction", 
#'                                        sample_size = 30, CI = 0.95)
#' }
#' @noRd
perform_unified_resampling = function(data_group, n_resamples, stat_fn, 
                                     original_input_is_vector, bootstrap_rarefaction,
                                     sample_size = NULL, CI = 0.95, ci_method = "percentile",
                                     parallel = FALSE, ncores = 1) {
  
  n_g = if (original_input_is_vector) length(data_group) else nrow(data_group)
  
  if (bootstrap_rarefaction == "bootstrap") {
    # Always use internal bootstrap to avoid Windows parallel issues
    observed_stat = stat_fn(data_group)
    resampled_values = perform_internal_bootstrap(data_group, n_resamples, stat_fn, 
                                                 original_input_is_vector, sample_size,
                                                 parallel, ncores)
    
    # Calculate confidence intervals
    if (CI == 1) {
      ci_result = c(min(resampled_values), max(resampled_values))
    } else if (requireNamespace("boot", quietly = TRUE) && ci_method == "bca") {
      # Create boot object for BCa intervals without using boot() function's parallel features
      ci_result = calculate_bca_ci(data_group, resampled_values, observed_stat, stat_fn, 
                                  original_input_is_vector, CI)
    } else {
      # Use percentile CI
      ci_result = calculate_percentile_ci(resampled_values, CI)
    }
    
  } else if (bootstrap_rarefaction == "rarefaction") {
    # Rarefaction without replacement
    if (is.null(sample_size)) {
      stop("sample_size must be provided for rarefaction")
    }
    
    resampled_values = perform_internal_rarefaction(data_group, n_resamples, sample_size,
                                                   stat_fn, original_input_is_vector,
                                                   parallel, ncores)
    observed_stat = mean(resampled_values)  # For rarefaction, observed is mean of resampled
    ci_result = calculate_percentile_ci(resampled_values, CI)
    
  } else if (bootstrap_rarefaction == "rarefaction_with_replacement") {
    # Rarefaction with replacement (essentially bootstrap with fixed sample size)
    if (is.null(sample_size)) {
      stop("sample_size must be provided for rarefaction_with_replacement")
    }
    
    observed_stat = stat_fn(data_group)
    resampled_values = perform_internal_bootstrap(data_group, n_resamples, stat_fn, 
                                                 original_input_is_vector, sample_size,
                                                 parallel, ncores)
    ci_result = calculate_percentile_ci(resampled_values, CI)
    
  } else {
    stop("bootstrap_rarefaction must be 'bootstrap', 'rarefaction', or 'rarefaction_with_replacement'")
  }
  
  return(list(
    resampled_values = resampled_values,
    observed_stat = observed_stat,
    CI_min = ci_result[1],
    CI_max = ci_result[2]
  ))
}


#' Internal bootstrap resampling function
#'
#' @param data_group Vector or matrix of data
#' @param n_resamples Number of resamples
#' @param stat_fn Statistic function
#' @param original_input_is_vector Logical for vector input
#' @param sample_size Optional sample size for resampling
#' @param parallel Logical for parallel processing
#' @param ncores Number of cores
#' @return Numeric vector of resampled values
#' @noRd
perform_internal_bootstrap = function(data_group, n_resamples, stat_fn, 
                                     original_input_is_vector, sample_size = NULL,
                                     parallel = FALSE, ncores = 1) {
  n_g = if (original_input_is_vector) length(data_group) else nrow(data_group)
  actual_sample_size = if (is.null(sample_size)) n_g else sample_size
  
  resample_fn = function(i) {
    sampled_idx = sample(seq_len(n_g), actual_sample_size, replace = TRUE)
    if (original_input_is_vector) {
      stat_fn(data_group[sampled_idx])
    } else {
      stat_fn(data_group[sampled_idx, , drop = FALSE])
    }
  }
  
  if (parallel && ncores > 1) {
    res_vals = unlist(safe_parallel_lapply(seq_len(n_resamples), resample_fn, ncores = ncores))
  } else {
    res_vals = unlist(lapply(seq_len(n_resamples), resample_fn))
  }
  
  return(res_vals)
}


#' Internal rarefaction resampling function
#'
#' @param data_group Vector or matrix of data
#' @param n_resamples Number of resamples
#' @param sample_size Sample size for rarefaction
#' @param stat_fn Statistic function
#' @param original_input_is_vector Logical for vector input
#' @param parallel Logical for parallel processing
#' @param ncores Number of cores
#' @return Numeric vector of resampled values
#' @noRd
perform_internal_rarefaction = function(data_group, n_resamples, sample_size, 
                                       stat_fn, original_input_is_vector,
                                       parallel = FALSE, ncores = 1) {
  n_g = if (original_input_is_vector) length(data_group) else nrow(data_group)
  
  if (sample_size > n_g) {
    stop("sample_size (", sample_size, ") cannot be larger than group size (", n_g, ")")
  }
  
  resample_fn = function(i) {
    sampled_idx = sample(seq_len(n_g), sample_size, replace = FALSE)
    if (original_input_is_vector) {
      stat_fn(data_group[sampled_idx])
    } else {
      stat_fn(data_group[sampled_idx, , drop = FALSE])
    }
  }
  
  if (parallel && ncores > 1) {
    res_vals = unlist(safe_parallel_lapply(seq_len(n_resamples), resample_fn, ncores = ncores))
  } else {
    res_vals = unlist(lapply(seq_len(n_resamples), resample_fn))
  }
  
  return(res_vals)
}


#' Calculate percentile confidence intervals
#'
#' @param resampled_values Numeric vector of resampled values
#' @param CI Confidence interval level
#' @return Named numeric vector with CI bounds
#' @noRd
calculate_percentile_ci = function(resampled_values, CI) {
  if (CI == 1) {
    return(c(min(resampled_values), max(resampled_values)))
  }
  
  alpha = (1 - CI) / 2
  ci_min = quantile(resampled_values, probs = alpha, names = FALSE, type = 7)
  ci_max = quantile(resampled_values, probs = 1 - alpha, names = FALSE, type = 7)
  
  return(c(ci_min, ci_max))
}


#' Calculate bias-corrected accelerated (BCa) confidence intervals
#'
#' Creates a boot object from existing bootstrap samples and uses boot.ci() 
#' to compute BCa intervals. This avoids Windows parallel issues with boot().
#'
#' @param data_group Original data (vector or matrix)
#' @param resampled_values Numeric vector of bootstrap replicates
#' @param observed_stat Observed statistic value
#' @param stat_fn Function used to compute statistic
#' @param original_input_is_vector Logical indicating vector input
#' @param CI Confidence interval level
#' @return Numeric vector with CI bounds
#' @noRd
calculate_bca_ci = function(data_group, resampled_values, observed_stat, stat_fn, 
                           original_input_is_vector, CI) {
  
  # Create a mock boot object compatible with boot.ci()
  n_resamples = length(resampled_values)
  n_g = if (original_input_is_vector) length(data_group) else nrow(data_group)
  
  # Generate the indices that would have been used (we don't have them from internal bootstrap)
  # For BCa calculation, boot.ci() needs the original indices, but we can approximate
  mock_indices = matrix(sample(seq_len(n_g), n_g * n_resamples, replace = TRUE), 
                       nrow = n_resamples, ncol = n_g)
  
  # Create boot object structure
  boot_obj = list(
    t0 = observed_stat,                    # Original statistic
    t = matrix(resampled_values, ncol = 1), # Bootstrap replicates
    R = n_resamples,                       # Number of replicates
    data = data_group,                     # Original data
    seed = .Random.seed,                   # For reproducibility
    statistic = function(data, indices) { # Statistic function
      if (original_input_is_vector) {
        stat_fn(data[indices])
      } else {
        stat_fn(data[indices, , drop = FALSE])
      }
    },
    sim = "ordinary",                      # Simulation type
    call = call("internal_bootstrap"),     # Call info
    stype = "i",                          # Index type
    strata = rep(1, n_g),                 # Strata info
    weights = rep(1/n_g, n_g)             # Weights
  )
  
  # Set class to make it a proper boot object
  class(boot_obj) = "boot"
  
  # Try to calculate BCa intervals
  tryCatch({
    ci_boot = boot::boot.ci(boot_obj, conf = CI, type = "bca")
    if (!is.null(ci_boot$bca)) {
      return(c(ci_boot$bca[4], ci_boot$bca[5]))  # BCa interval
    } else {
      # Fallback to percentile if BCa calculation fails
      warning("BCa calculation failed, using percentile intervals")
      return(calculate_percentile_ci(resampled_values, CI))
    }
  }, error = function(e) {
    # Fallback to percentile if any error occurs
    warning("BCa calculation failed with error: ", e$message, ". Using percentile intervals.")
    return(calculate_percentile_ci(resampled_values, CI))
  })
}

