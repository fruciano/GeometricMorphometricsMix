# NOTE: This file contains internal helper functions for fitting PGLS
# and generating simulations from fitted PGLS models. The package
# DESCRIPTION should list the following packages in the Imports field:
#   Imports:
#     phylolm
#     mvMORPH
# The functions below are internal (not exported) and are intended to be
# used by other functions in this package.

#' Fit a phylogenetic GLS model (univariate or multivariate)
#'
#' Internal helper that dispatches to phylolm::phylolm for univariate
#' responses and mvMORPH::mvgls for multivariate responses.
#'
#' @param Y Response: vector, matrix or data.frame (tips in rows).
#' @param X Predictor: numeric vector or matrix (length/rows match tree tips).
#' @param tree A phylo object with tip labels matching rows of Y/X.
#' @param model Character, evolutionary model to pass to the underlying fitter (default "BM").
#' @return A fitted model object (class 'phylolm' or 'mvgls').
#' @keywords internal
#' @noRd
fit_pgls = function(tree, Y, X, model = "BM"){

  # Decide whether Y is univariate or multivariate and call the
  # appropriate fitting function using namespace-qualified calls.
  is_univariate = is.vector(Y) || (is.matrix(Y) && ncol(Y) == 1) || (is.data.frame(Y) && ncol(Y) == 1)

  if(is_univariate){
    # coerce to numeric vector and ensure rownames
    if(is.matrix(Y) || is.data.frame(Y)){
      yvec = as.numeric(Y[,1])
      rn = rownames(Y)
    } else {
      yvec = as.numeric(Y)
      rn = names(Y)
    }
    if(is.null(rn)) rn = tree$tip.label
    names(yvec) = rn
    dat = data.frame(Y = yvec, X = X)
    rownames(dat) = rn
    
    # call phylolm for univariate response
  PGLS_model = phylolm::phylolm(formula = Y ~ X, data = dat, phy = tree, model = ifelse(missing(model), "BM", model))
    PGLS_model$tree = tree  # store tree in the model for later use
  } else {
    # multivariate: ensure matrix and rownames
    Ymat = if(is.data.frame(Y)) as.matrix(Y) else Y
    if(is.null(rownames(Ymat))) rownames(Ymat) = tree$tip.label

    # call mvgls for multivariate response
    PGLS_model = mvMORPH::mvgls(formula = Ymat ~ X, tree = tree, model = model)
  }

  return(PGLS_model)
}


#' Generate simulated datasets from a fitted PGLS model
#'
#' Generic wrapper that routes to the appropriate simulator depending on
#' whether the fitted model is univariate (phylolm) or multivariate (mvgls).
#'
#' @param PGLS_model Fitted model from fit_pgls / phylolm::phylolm or mvMORPH::mvgls.
#' @param nsim Number of simulated datasets to produce.
#' @param model Evolutionary model name passed to mvMORPH::mvSIM (e.g. "BM").
#' @return A list of simulated response matrices/vectors (length = nsim).
#' @keywords internal
#' @noRd
PGLS_sim_gen = function(PGLS_model, nsim = 1000, model = "BM"){
  if(inherits(PGLS_model, "phylolm")){
    # Univariate case
    univariate_PGLS_sim_gen(PGLS_model, nsim, model)
  } else if (inherits(PGLS_model, "mvgls")){
    # Multivariate case: map BM -> BM1 for mvMORPH when needed
    if(model == "BM") model = "BM1"
    multivariate_PGLS_sim_gen(PGLS_model, nsim, model)
  } else {
    stop("PGLS_model must be of class 'phylolm' or 'mvgls'")
  }
}


#' Simulate multivariate datasets from a fitted mvgls/phylolm object
#'
#' Uses the residual covariance from the fitted model as the sigma parameter
#' for mvMORPH::mvSIM, and scales simulated errors so that their total
#' residual variance matches the fitted model.
#'
#' @param PGLS_model A fitted mvgls (or compatible) model with $residuals and a stored tree.
#' @param nsim Number of simulations to perform.
#' @param model Model name for mvMORPH::mvSIM (e.g. "BM1").
#' @return List of simulated response matrices (each matrix has rows = tips, cols = traits).
#' @keywords internal
#' @noRd
multivariate_PGLS_sim_gen = function(PGLS_model, nsim = 1000, model = "BM1"){
  if(missing(PGLS_model) || is.null(PGLS_model$residuals)){
    stop("PGLS_model must be provided and contain residuals")
  }

  # try to extract tree from common locations
  tree = NULL
  if(!is.null(PGLS_model$variables) && !is.null(PGLS_model$variables$tree)){
    tree = PGLS_model$variables$tree
  }
  if(is.null(tree) && !is.null(PGLS_model$phy)) tree = PGLS_model$phy
  if(is.null(tree) && !is.null(PGLS_model$tree)) tree = PGLS_model$tree
  if(is.null(tree)) stop("Could not find tree in PGLS_model (looked for $variables$tree, $phy, or $tree)")

  # covariance for simulation: use residual covariance from the fitted model
  cov_mat = stats::cov(PGLS_model$residuals)

  Sims_error = suppressMessages(suppressWarnings(
    suppressWarnings(capture.output(
      mvMORPH::mvSIM(tree, nsim = nsim, model = model, param = list(sigma = cov_mat)), file = NULL))))
  # Compute simulations using mvMORPH::mvSIM()


  # traces and cross-trace (use fitted and residuals if available)
  traceR = sum(diag(stats::cov(PGLS_model$residuals)))
  crossTrace = if(!is.null(PGLS_model$fitted) && !is.null(PGLS_model$residuals))
    2 * sum(diag(stats::cov(PGLS_model$fitted, PGLS_model$residuals))) else 0

  # required residual variance in sims
  target_tot_resid_var = traceR + crossTrace

  # now scale simulated errors to this target
  Sims_error = lapply(Sims_error, function(X){
    Xmat = as.matrix(X)
    sim_var = sum(diag(stats::cov(Xmat)))
    if(sim_var == 0) return(Xmat)
    scaling_factor = sqrt(target_tot_resid_var / sim_var)
    Xmat * scaling_factor
  })

  # add fitted values if available, otherwise return scaled errors
  if(!is.null(PGLS_model$fitted)){
    Model_plus_error = lapply(Sims_error, function(X) PGLS_model$fitted + X)
  } else {
    Model_plus_error = Sims_error
  }

  return(Model_plus_error)
}


#' Simulate univariate datasets from a fitted phylolm object
#'
#' Reconstructs the observed response as fitted + residuals and simulates
#' residuals under the specified evolutionary model. Simulated residuals
#' are scaled so that total residual variance matches the fitted model.
#'
#' @param PGLS_model Fitted phylolm model (must contain $fitted.values and $residuals).
#' @param nsim Number of simulations to produce.
#' @param model Evolutionary model name for mvMORPH::mvSIM ("BM" maps to "BM1").
#' @return List of simulated response vectors (each numeric vector of length = number of tips).
#' @keywords internal
#' @noRd
univariate_PGLS_sim_gen = function(PGLS_model, nsim = 1000, model = "BM"){
  # Extract fitted values and residuals
  Fitted = as.numeric(PGLS_model$fitted.values)
  R_emp = as.numeric(PGLS_model$residuals)

  # Reconstruct observed response as fitted + residuals
  Yvec = Fitted + R_emp

  # compute sample variances and covariance (R uses n-1 denom)
  varY   = stats::var(Yvec)
  varF   = stats::var(Fitted)
  varR   = stats::var(R_emp)
  covFR  = stats::cov(Fitted, R_emp)

  tree = NULL
  if(!is.null(PGLS_model$tree)) tree = PGLS_model$tree
  if(is.null(tree) && !is.null(PGLS_model$phy)) tree = PGLS_model$phy
  if(is.null(tree) && !is.null(PGLS_model$variables) && !is.null(PGLS_model$variables$tree)) tree = PGLS_model$variables$tree
  if(is.null(tree)) stop("Could not find tree in PGLS_model for simulation")

  # two equivalent ways to get target residual variance
  target_var1 = varR + 2 * covFR
  target_var2 = varY - varF

  # choose target (prefer varY - varF)
  target_resid_var = target_var2

  # safeguard: if negative or non-finite, fallback to matching varR
  if(!is.finite(target_resid_var) || target_resid_var <= 0){
    warning("Computed target residual variance not positive; falling back to empirical residual variance.")
    target_resid_var = varR
  }

  if (model == "BM") model = "BM1"

  Sims_error = suppressMessages(suppressWarnings(
    suppressWarnings(capture.output(
      mvMORPH::mvSIM(tree, nsim = nsim, model = model, param = list(sigma = target_resid_var)), file = NULL))))

  # mvSIM for univariate may return a matrix with columns = simulations
  if(is.matrix(Sims_error)){
    Sims_error_scaled = apply(Sims_error, 2, function(x){
      sim_var = stats::var(as.numeric(x))
      if(sim_var <= 0 || !is.finite(sim_var)){
        stop("Simulated residual variance is zero or not finite.")
      }
      scaling = sqrt(target_resid_var / sim_var)
      as.numeric(x) * scaling
    })

    # create simulated datasets
    Model_plus_error = lapply(seq(ncol(Sims_error_scaled)), function(i) as.numeric(Fitted + Sims_error_scaled[,i]))
  } else if(is.list(Sims_error)){
    # If mvSIM returned a list of vectors/matrices
    Model_plus_error = lapply(Sims_error, function(x){
      xvec = as.numeric(x)
      sim_var = stats::var(xvec)
      if(sim_var <= 0 || !is.finite(sim_var)) stop("Simulated residual variance is zero or not finite.")
      scaling = sqrt(target_resid_var / sim_var)
      as.numeric(Fitted + xvec * scaling)
    })
  } else {
    stop("Unexpected return type from mvMORPH::mvSIM for univariate simulation")
  }

  return(Model_plus_error)
}

