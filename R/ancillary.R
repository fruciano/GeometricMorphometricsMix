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
#' This function is used internally by a few functions in this
#' package but is exported so that users can reuse it for lightweight
#' cross-platform parallel work (for example on Windows where forking
#' is not available).
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
#' # more complex example: simulate work with a pause to compare
#' # serial vs parallel timing
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
safe_parallel_lapply = function(
  X, FUN, ncores = (parallel::detectCores(logical = FALSE)-1),
  packages = NULL, export = NULL, seed = NULL,
  type = c("auto", "mclapply", "psock")
) {
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

    if (!is.null(export)) parallel::clusterExport(cl, export,
                                                   envir = parent.frame())
    if (!is.null(packages)) {
      # export the packages vector so workers can access it
      parallel::clusterExport(cl, varlist = "packages",
                              envir = parent.frame())
      parallel::clusterEvalQ(cl, {
        for (pk in packages) {
          try(loadNamespace(pk), silent = TRUE)
        }
        NULL
      })
    }
    if (!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)

    # run in a tryCatch to ensure cluster is stopped and errors are
    # informative
    res = tryCatch({
      parallel::parLapply(cl, X, FUN)
    }, error = function(e) {
      stop(
        "Error while running parLapply on PSOCK cluster: ",
        conditionMessage(e)
      )
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
    res = parallel::mclapply(
      X, FUN, mc.cores = ncores, mc.preschedule = TRUE,
      mc.set.seed = TRUE
    )
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
      parallel::mclapply(
        X, FUN, mc.cores = ncores, mc.preschedule = TRUE,
        mc.set.seed = TRUE
      )
    }, silent = TRUE)
    if (!inherits(try_res, "try-error")) return(try_res)
    # fallback
    return(run_psock())
  }
}

