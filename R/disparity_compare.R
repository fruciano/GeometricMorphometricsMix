#' Compare disparity between two groups
#'
#' Compares disparity using permutation, within-group bootstrap, or
#' simulations from a fitted phylogenetic model.
#'
#' @param Data Numeric matrix or data frame, with observations in rows and
#'   variables in columns. Measurements must be finite. Phylogenetic analyses
#'   require at least two variables and unique species names as row names.
#' @param group Factor, or vector coercible to a factor, with one entry per
#'   row of \code{Data}. Exactly two nonmissing groups must be represented,
#'   each with at least two observations. Unused levels are dropped.
#' @param statistic Character vector specifying one or both of
#'   \code{"multivariate_variance"} (sum of sample variances) and
#'   \code{"mean_pairwise_euclidean_distance"}.
#' @param comparison One of \code{"log_ratio"} (default), \code{"ratio"},
#'   \code{"difference"}, or \code{"absolute_difference"}.
#' @param method One of \code{"permutation"} (default), \code{"bootstrap"},
#'   or \code{"phylogenetic"}.
#' @param n_resamples Positive integer number of replicates (default 999).
#'   Bootstrap and phylogenetic methods require at least two.
#' @param CI Central interval level, strictly between zero and one, for
#'   bootstrap confidence intervals or simulated/permuted null intervals.
#' @param tree A \code{"phylo"} tree, required only for the phylogenetic
#'   method. Its unique tip labels must match the full \code{Data} row names
#'   exactly. No pruning is performed.
#' @param ... Named arguments passed to \code{mvMORPH::mvgls()}, such as
#'   \code{model="OU"} or \code{REML=FALSE}. The default model is \code{"BM"};
#'   fitting uses the default method of \code{mvgls()}. Formula, response,
#'   data, and tree are controlled internally. A non-NULL \code{tree} and
#'   arguments in \code{...} are rejected for nonphylogenetic methods.
#'
#' @details
#' Both measures describe within-group disparity and are translation
#' invariant. No measurement adjustment or Procrustes fit is performed.
#' All requested statistics share each replicate.
#'
#' @section Groups and missing membership:
#' Factor level order is retained; other vectors use \code{factor()} order.
#' Differences are group 1 minus group 2; ratios are group 1 divided by
#' group 2, with group 1 the first level.
#'
#' \code{NA} membership excludes an observation from disparity calculations
#' and, for permutation/bootstrap, sampling. For phylogenetic analyses it
#' remains in fitting and simulation: the focal groups may be subsets of
#' the full dataset. Missing measurements are not supported.
#'
#' @section Resampling methods:
#' \describe{
#'   \item{Permutation}{Pools focal observations and shuffles them while
#'     preserving group sizes.}
#'   \item{Bootstrap}{Samples independently with replacement within each
#'     group, preserving its original size.}
#'   \item{Phylogenetic}{Fits an intercept-only multivariate GLS model to all
#'     observations with \code{mvMORPH::mvgls()}, then simulates using
#'     \code{stats::simulate()}. Group membership and phylogenetic positions
#'     remain fixed. Inference is conditional on the fitted model
#'     (Brownian motion as default)and its
#'     estimated parameters; fitting uncertainty is not propagated.}
#' }
#'
#' @section Statistical inference:
#' Permutation and phylogenetic analyses use the upper empirical tail for
#' absolute differences; other comparisons use twice the smaller tail,
#' capped at one. Ties are included, with the plus-one correction in both
#' tails. P-values are unadjusted for multiple statistics.
#'
#' Bootstrap produces no p-value: percentile confidence intervals describe
#' sampling uncertainty. Permutation/phylogenetic intervals describe null
#' distributions, which need not be centered on equal disparity.
#'
#' @note
#' \itemize{
#'   \item Permutation requires exchangeability. Pooling groups with different
#'     locations can affect the null distribution despite translation-invariant
#'     statistics. Consider subtracting each original group's column means
#'     when comparing dispersion in groups with different means. Centering
#'     is not automatic and does not guarantee calibration under all
#'     distributional or sample-size differences.
#'   \item Group translations leave bootstrap contrasts unchanged. Centering
#'     is not recommended for phylogenetic analyses: it changes the fitting
#'     data and potentially the simulated null distribution.
#'   \item Ratios require a positive denominator; log ratios require both
#'     disparities positive and use the difference of logs. Undefined observed
#'     contrasts cause an error. Undefined replicate contrasts remain NA;
#'     any occurrence triggers a warning and withholds all inferential
#'     summaries for that statistic. Observed values and valid counts remain
#'     available; replicates are neither replaced nor silently discarded.
#' }
#'
#' @return A list of class \code{"disparity_compare"}, displayed compactly by
#'   \code{\link{print.disparity_compare}}, with elements:
#' \describe{
#'   \item{results}{Data frame with one row per statistic:
#'     \describe{
#'       \item{Common columns}{\code{statistic}, \code{Observed_grp1},
#'         \code{Observed_grp2}, \code{contrast}, \code{p_value}, and
#'         \code{n_valid}. Bootstrap p-values are NA.}
#'       \item{Bootstrap columns}{\code{bootstrap_mean}, \code{bootstrap_sd},
#'         \code{bootstrap_bias}, \code{CI_lower}, \code{CI_median},
#'         \code{CI_upper}, \code{null_value}, and
#'         \code{interval_contains_null}. Inclusion is NA for absolute
#'         differences, where inclusion of zero is uninformative.}
#'       \item{Null-model columns}{\code{null_mean}, \code{null_sd},
#'         \code{null_lower}, \code{null_median}, \code{null_upper}, and
#'         \code{SES}: observed contrast minus null mean, divided by null SD.
#'         SES is NA when that SD is zero or cannot be estimated.}
#'     }}
#'   \item{replicate_contrasts}{Numeric matrix: replicates in rows,
#'     statistics in columns, including NA for undefined contrasts.}
#'   \item{group_levels, group_sizes}{Labels in comparison order and sample sizes.}
#'   \item{statistic, comparison, method, n_resamples, CI_level}{Settings used.}
#'   \item{evolutionary_model, fit}{Model name and fitted \code{mvgls} object
#'     for phylogenetic analyses; otherwise NULL.}
#'   \item{call}{The matched function call.}
#' }
#'
#' @seealso \code{\link{disparity_test}}, \code{\link{disparity_resample}}
#' @references Foote M. 1997. The evolution of morphological diversity.
#'   Annual Review of Ecology and Systematics 28:129-152.
#' @references Efron B, Tibshirani RJ. 1993. An Introduction to the Bootstrap.
#'   Chapman & Hall.
#' @references Clavel J, Aristide L, Morlon H. 2019. A penalized likelihood
#'   framework for high-dimensional phylogenetic comparative methods and
#'   an application to New-World monkeys brain evolution.
#'   Systematic Biology 68:93-116.
#'
#' @examples
#' set.seed(123)
#' Data=matrix(rnorm(160), nrow=40)
#' grp=factor(rep(c("A", "B"), each=20), levels=c("A", "B"))
#' disparity_compare(Data, grp, n_resamples=99)
#'
#' # Request both measures, sharing the same permutations
#' disparity_compare(
#'   Data, grp, statistic=c("multivariate_variance",
#'                         "mean_pairwise_euclidean_distance"),
#'   comparison="absolute_difference", n_resamples=99
#' )$results
#'
#' disparity_compare(Data, grp, method="bootstrap", n_resamples=199)
#'
#' # Nonfocal species remain in the phylogenetic fit
#' \dontrun{
#' if (requireNamespace("ape", quietly=TRUE) &&
#'     requireNamespace("mvMORPH", quietly=TRUE)) {
#'   tree=ape::rtree(40)
#'   rownames(Data)=tree$tip.label
#'   grp_phylo=grp
#'   grp_phylo[c(1:5, 21:25)]=NA
#'   disparity_compare(Data, grp_phylo, method="phylogenetic",
#'                     tree=tree, n_resamples=199, model="BM")$results
#' }
#' }
#' @export
disparity_compare=function(Data, group,
                           statistic="multivariate_variance",
                           comparison="log_ratio", method="permutation",
                           n_resamples=999, CI=0.95, tree=NULL, ...) {
  function_call=match.call()
  # Validate data and the two focal groups
  if (!(is.matrix(Data) || is.data.frame(Data))) {
    stop("'Data' must be a numeric matrix or data frame.")
  }
  Data=as.matrix(Data)
  if (!is.numeric(Data) || ncol(Data) < 1 || any(!is.finite(Data))) {
    stop("'Data' must have at least one column of finite numeric measurements.")
  }
  if (!is.atomic(group) || !is.null(dim(group)) ||
      length(group) != nrow(Data)) {
    stop("'group' must be a vector with one entry per row of 'Data'.")
  }
  group=factor(group, exclude=NA)
  group_levels=levels(group)
  if (length(group_levels) != 2) {
    stop("'group' must contain exactly two observed nonmissing groups.")
  }
  group_sizes=table(group)
  if (any(group_sizes < 2)) {
    stop("Each group must contain at least two observations.")
  }
  group_sizes=setNames(as.integer(group_sizes), group_levels)

  allowed_statistics=c("multivariate_variance",
                       "mean_pairwise_euclidean_distance")
  if (!is.character(statistic) || length(statistic) < 1 ||
      anyNA(statistic) || any(!statistic %in% allowed_statistics) ||
      anyDuplicated(statistic)) {
    stop("'statistic' must specify one or both supported statistics, without duplicates.")
  }
  comparison=match.arg(comparison, c("log_ratio", "ratio", "difference",
                                    "absolute_difference"))
  method=match.arg(method, c("permutation", "bootstrap", "phylogenetic"))
  if (!is.numeric(n_resamples) || length(n_resamples) != 1 ||
      !is.finite(n_resamples) || n_resamples < 1 ||
      n_resamples > .Machine$integer.max || n_resamples != floor(n_resamples)) {
    stop("'n_resamples' must be a positive integer.")
  }
  n_resamples=as.integer(n_resamples)
  if (method != "permutation" && n_resamples < 2) {
    stop("Bootstrap and phylogenetic methods require at least two replicates.")
  }
  if (!is.numeric(CI) || length(CI) != 1 || !is.finite(CI) ||
      CI <= 0 || CI >= 1) {
    stop("'CI' must be strictly between zero and one.")
  }
  dots=list(...)
  if (method != "phylogenetic" && (!is.null(tree) || length(dots) > 0)) {
    stop("'tree' and arguments in ... are only used for the phylogenetic method.")
  }

  # Keep the full dataset for the phylogenetic model
  if (method == "phylogenetic") {
    if (!requireNamespace("mvMORPH", quietly=TRUE)) {
      stop("Package 'mvMORPH' is required for the phylogenetic method.")
    }
    if (!inherits(tree, "phylo")) {
      stop("'tree' must be an object of class 'phylo'.")
    }
    if (ncol(Data) < 2) {
      stop("The phylogenetic method requires at least two variables.")
    }
    species=rownames(Data)
    tips=tree$tip.label
    if (is.null(species) || anyNA(species) || any(species == "") ||
        anyDuplicated(species) || is.null(tips) || anyNA(tips) ||
        any(tips == "") || anyDuplicated(tips) || !setequal(species, tips)) {
      stop("Unique, nonempty row names of 'Data' must match tree tip labels exactly.")
    }
    if (length(dots) > 0 &&
        (is.null(names(dots)) || any(names(dots) == "") ||
         anyDuplicated(names(dots)))) {
      stop("Arguments in ... must have unique names.")
    }
    if (any(names(dots) %in% c("formula", "data", "tree", "response"))) {
      stop("The formula, data, tree and response are controlled internally.")
    }
    if (!"model" %in% names(dots)) dots$model="BM"
    order_idx=match(tips, species)
    Data=Data[order_idx,,drop=FALSE]
    group=group[order_idx]
  } else {
    keep=which(!is.na(group))
    Data=Data[keep,,drop=FALSE]
    group=group[keep]
  }
  idx1=which(group == group_levels[1])
  idx2=which(group == group_levels[2])
  X1=Data[idx1,,drop=FALSE]
  X2=Data[idx2,,drop=FALSE]
  n1=nrow(X1)
  n2=nrow(X2)

  # Small local functions share the computations across all methods
  compute_disparity=function(X) {
    vapply(statistic, function(s) {
      if (s == "multivariate_variance") muvar(X) else meanpairwiseEuclideanD(X)
    }, numeric(1))
  }
  compute_contrast=function(d1, d2) {
    contrast=rep(NA_real_, length(statistic))
    valid=is.finite(d1) & is.finite(d2)
    if (comparison == "difference") contrast[valid]=d1[valid]-d2[valid]
    if (comparison == "absolute_difference") {
      contrast[valid]=abs(d1[valid]-d2[valid])
    }
    if (comparison == "ratio") {
      valid=valid & d2 > 0
      contrast[valid]=d1[valid]/d2[valid]
    }
    if (comparison == "log_ratio") {
      valid=valid & d1 > 0 & d2 > 0
      contrast[valid]=log(d1[valid])-log(d2[valid])
    }
    contrast[!is.finite(contrast)]=NA_real_
    return(contrast)
  }
  observed1=compute_disparity(X1)
  observed2=compute_disparity(X2)
  observed=compute_contrast(observed1, observed2)
  if (anyNA(observed)) {
    stop("Undefined observed contrast for: ",
         paste(statistic[is.na(observed)], collapse=", "),
         ". Ratios need a positive denominator; log ratios need positive disparities.")
  }

  fit=NULL
  evolutionary_model=NULL
  simulations=NULL
  if (method == "phylogenetic") {
    fit=do.call(mvMORPH::mvgls,
                c(list(formula=Y ~ 1, data=list(Y=Data), tree=tree), dots))
    evolutionary_model=fit$model
    simulations=stats::simulate(fit, nsim=n_resamples)
    if (!is.list(simulations) || length(simulations) != n_resamples) {
      stop("Expected a list of simulated datasets from simulate.mvgls().")
    }
  }

  # Compute contrasts without retaining permuted or bootstrapped datasets
  replicate_contrasts=matrix(NA_real_, nrow=n_resamples,
                             ncol=length(statistic),
                             dimnames=list(NULL, statistic))
  for (i in seq_len(n_resamples)) {
    if (method == "permutation") {
      # Preserve the full-permutation draw used by disparity_test()
      perm_idx=sample(seq_len(nrow(Data)), nrow(Data), replace=FALSE)
      draw1=Data[perm_idx[seq_len(n1)],,drop=FALSE]
      draw2=Data[perm_idx[(n1+1):(n1+n2)],,drop=FALSE]
    } else if (method == "bootstrap") {
      draw1=X1[sample(seq_len(n1), n1, replace=TRUE),,drop=FALSE]
      draw2=X2[sample(seq_len(n2), n2, replace=TRUE),,drop=FALSE]
    } else {
      simulated=as.matrix(simulations[[i]])
      if (!is.numeric(simulated) || !identical(dim(simulated), dim(Data))) {
        stop("Unexpected dimensions or nonnumeric data in simulated dataset ", i, ".")
      }
      # mvgls simulations follow the model's species order when unnamed
      if (is.null(rownames(simulated))) rownames(simulated)=tree$tip.label
      if (anyDuplicated(rownames(simulated)) ||
          !setequal(rownames(simulated), tree$tip.label)) {
        stop("Species labels in simulated dataset ", i, " do not match the tree.")
      }
      simulated=simulated[tree$tip.label,,drop=FALSE]
      draw1=simulated[idx1,,drop=FALSE]
      draw2=simulated[idx2,,drop=FALSE]
    }
    replicate_contrasts[i,]=compute_contrast(compute_disparity(draw1),
                                             compute_disparity(draw2))
  }

  # Keep confidence intervals distinct from intervals under a null model
  labels=c(multivariate_variance="Multivariate variance",
           mean_pairwise_euclidean_distance="Mean pairwise Euclidean distance")
  Results=data.frame(statistic=statistic,
                     Observed_grp1=unname(observed1),
                     Observed_grp2=unname(observed2), contrast=observed,
                     p_value=NA_real_,
                     n_valid=colSums(is.finite(replicate_contrasts)),
                     row.names=unname(labels[statistic]))
  if (method == "bootstrap") {
    Results$bootstrap_mean=NA_real_
    Results$bootstrap_sd=NA_real_
    Results$bootstrap_bias=NA_real_
    Results$CI_lower=NA_real_
    Results$CI_median=NA_real_
    Results$CI_upper=NA_real_
    Results$null_value=if (comparison == "ratio") 1 else 0
    Results$interval_contains_null=NA
  } else {
    Results$null_mean=NA_real_
    Results$null_sd=NA_real_
    Results$null_lower=NA_real_
    Results$null_median=NA_real_
    Results$null_upper=NA_real_
    Results$SES=NA_real_
  }
  probs=c((1-CI)/2, 0.5, 1-(1-CI)/2)
  for (j in seq_along(statistic)) {
    if (Results$n_valid[j] != n_resamples) {
      warning("Undefined replicate contrasts for ", statistic[j], " (",
              Results$n_valid[j], " of ", n_resamples,
              " valid). Inferential summaries are withheld.", call.=FALSE)
      next
    }
    values=replicate_contrasts[,j]
    interval=stats::quantile(values, probs=probs, names=FALSE)
    if (method == "bootstrap") {
      Results$bootstrap_mean[j]=mean(values)
      Results$bootstrap_sd[j]=stats::sd(values)
      Results$bootstrap_bias[j]=mean(values)-observed[j]
      Results$CI_lower[j]=interval[1]
      Results$CI_median[j]=interval[2]
      Results$CI_upper[j]=interval[3]
      if (comparison != "absolute_difference") {
        Results$interval_contains_null[j]=interval[1] <= Results$null_value[j] &&
          interval[3] >= Results$null_value[j]
      }
    } else {
      p_upper=(1+sum(values >= observed[j]))/(n_resamples+1)
      p_lower=(1+sum(values <= observed[j]))/(n_resamples+1)
      Results$p_value[j]=if (comparison == "absolute_difference") {
        p_upper
      } else min(1, 2*min(p_lower, p_upper))
      Results$null_mean[j]=mean(values)
      Results$null_sd[j]=stats::sd(values)
      Results$null_lower[j]=interval[1]
      Results$null_median[j]=interval[2]
      Results$null_upper[j]=interval[3]
      if (is.finite(Results$null_sd[j]) && Results$null_sd[j] > 0) {
        Results$SES[j]=(observed[j]-mean(values))/Results$null_sd[j]
      }
    }
  }
  result=list(results=Results, replicate_contrasts=replicate_contrasts,
              group_levels=group_levels, group_sizes=group_sizes,
              statistic=statistic, comparison=comparison, method=method,
              n_resamples=n_resamples, CI_level=CI,
              evolutionary_model=evolutionary_model, fit=fit,
              call=function_call)
  class(result)=c("disparity_compare", "list")
  return(result)
}

#' Print the results of a disparity comparison
#'
#' @param x An object returned by \code{\link{disparity_compare}}.
#' @param digits Number of significant digits used for display.
#' @param ... Additional arguments passed to data-frame printing.
#' @return The original object, invisibly.
#' @export
print.disparity_compare=function(x, digits=4, ...) {
  method_label=switch(x$method,
                      permutation="Permutation",
                      bootstrap="Within-group bootstrap",
                      phylogenetic="Phylogenetic simulation")
  direction=switch(x$comparison,
                   difference="group 1 - group 2",
                   absolute_difference="abs(group 1 - group 2)",
                   ratio="group 1 / group 2",
                   log_ratio="ln(group 1 / group 2)")
  cat("Disparity comparison\n")
  cat("Method: ", method_label, "\n", sep="")
  cat("Group 1: ", x$group_levels[1], " (n=", x$group_sizes[1], ")\n", sep="")
  cat("Group 2: ", x$group_levels[2], " (n=", x$group_sizes[2], ")\n", sep="")
  cat("Comparison: ", direction, "\n", sep="")
  cat("Requested replicates: ", x$n_resamples, "\n", sep="")
  if (x$method == "phylogenetic") {
    cat("Evolutionary model: ", x$evolutionary_model, "\n", sep="")
  }

  observed=x$results[,c("Observed_grp1", "Observed_grp2", "contrast"),drop=FALSE]
  colnames(observed)=c("Group 1", "Group 2", "Contrast")
  cat("\nObserved disparities and contrast:\n")
  print(observed, digits=digits, ...)

  if (x$method == "bootstrap") {
    inference=x$results[,c("CI_lower", "CI_upper"),drop=FALSE]
    colnames(inference)=c("Lower", "Upper")
    if (x$comparison != "absolute_difference") {
      inference[["Includes equality"]]=x$results$interval_contains_null
    }
    cat("\n", format(100*x$CI_level, digits=digits, trim=TRUE),
        "% bootstrap confidence intervals:\n", sep="")
  } else {
    inference=data.frame("p-value"=format.pval(x$results$p_value, digits=digits),
                          row.names=rownames(x$results), check.names=FALSE)
    cat("\nInference:\n")
  }
  inference[["Valid/Requested"]]=paste0(x$results$n_valid, "/", x$n_resamples)
  print(inference, digits=digits, ...)
  invisible(x)
}
