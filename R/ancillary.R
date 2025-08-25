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
#' A small helper that selects an appropriate parallel backend depending on the
#' platform. On Windows it uses a PSOCK cluster (parLapply), on Unix-like
#' systems it tries to use forking via parallel::mclapply and falls back to a
#' PSOCK cluster if forking is not available or fails. If \code{workers==1}
#' this simply calls \code{lapply}. The function supports exporting objects and
#' loading packages on worker processes and has basic RNG handling via
#' \code{seed}.
#'
#' This helper is intended for internal package use (not exported).
#'
#' @noRd
safe_parallel_lapply <- function(X, FUN, ncores = parallel::detectCores(logical = FALSE),
                                 packages = NULL, export = NULL, seed = NULL,
                                 type = c("auto", "mclapply", "psock")) {
  type <- match.arg(type)
  ncores <- as.integer(max(1, min(length(X), ncores)))

  if (ncores <= 1) return(lapply(X, FUN))

  os_type <- .Platform$OS.type  # "windows" or "unix"

  # Helper to run on PSOCK cluster
  run_psock <- function() {
  cl <- parallel::makeCluster(ncores, type = "PSOCK")
    # ensure cluster is stopped on exit of this function
    on.exit(parallel::stopCluster(cl), add = TRUE)

    if (!is.null(export)) parallel::clusterExport(cl, export, envir = parent.frame())
    if (!is.null(packages)) {
      parallel::clusterEvalQ(cl, {
        invisible(lapply(packages, function(pk) require(pk, character.only = TRUE)))
      })
    }
    if (!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)

    res <- parallel::parLapply(cl, X, FUN)
    return(res)
  }

  # If user forced a backend
  if (type == "mclapply") {
    if (os_type == "windows") {
      # cannot fork on windows; fallback to psock
      return(run_psock())
    }
    if (!is.null(seed)) set.seed(seed)
    return(parallel::mclapply(X, FUN, mc.cores = ncores, mc.preschedule = TRUE, mc.set.seed = TRUE))
  }

  if (type == "psock") return(run_psock())

  # auto: choose best available
  if (os_type == "windows") {
    return(run_psock())
  } else {
    # try mclapply but fallback to PSOCK on error
    try_res <- try({
      if (!is.null(seed)) set.seed(seed)
      parallel::mclapply(X, FUN, mc.cores = ncores, mc.preschedule = TRUE, mc.set.seed = TRUE)
    }, silent = TRUE)
    if (!inherits(try_res, "try-error")) return(try_res)
    # fallback
    return(run_psock())
  }
}

