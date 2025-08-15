#' Phylogenetic Weighted Regression for Comparative Data
#'
#' Performs a phylogenetic weighted regression for each tip in a phylogenetic tree,
#' using a Brownian motion correlation structure. The function fits a weighted linear model
#' for each tip, computes R-squared, intercepts, and coefficients, and returns these along
#' with the weights matrix and the input tree.
#'
#' @param tree An object of class `phylo` (from the `ape` package) representing the phylogenetic tree.
#' @param Y A numeric vector, matrix, or data frame of response variables. Each column is treated as a trait.
#' @param X A numeric vector or factor of predictor values (e.g., trait or grouping variable).
#'
#' @return A list containing:
#'   \item{tree}{The input phylogenetic tree.}
#'   \item{rsquared}{A named vector of R-squared values for each tip.}
#'   \item{intercepts}{A matrix of intercepts for each tip.}
#'   \item{coefficients}{A matrix of regression coefficients for each tip.}
#'   \item{weights}{The phylogenetic correlation matrix used as weights.}
#'
#' @details
#' This function uses `ape::corBrownian` to compute the phylogenetic correlation structure,
#' and `nlme::corMatrix` to extract the weights matrix. Weighted linear models are fit using
#' `stats::lm`. The function is suitable for comparative analyses where phylogenetic signal
#' and trait evolution are modeled under Brownian motion.
#'
#' The returned object has class `phylo_wregr_fit` and has an associated S3 plot method
#' for visualizing results (see \code{\link{plot.phylo_wregr_fit}}).
#'
#' @importFrom ape corBrownian
#' @importFrom nlme corMatrix Initialize
#' @importFrom stats lm var prcomp
#' @importFrom mclust Mclust
#'
#' @note This function may give problems (NA in the output)
#' in the case of a single species sister to all the rest of the taxa
#' (this is due to the fact that all the phylogenetic correlations with
#' other taxa, on which the weighting is base, would be 0)
#'
#' @seealso \code{\link{plot.phylo_wregr_fit}} for plotting the results
#'
#' @examples
#' # Simulate a tree and traits
#' library(phytools)
#'
#' set.seed(123)
#'
#' tree = pbtree(n=20)
#' # Generate a tre with 20 tips
#'
#' Y = matrix(rnorm(100), nrow=20)
#' X = rnorm(20)
#' result = phylo_w_regr_fit(tree, Y, X)
#' str(result)
#'
#' # Plot R-squared values with phylogeny
#' plot(result, what = "rsquared")
#'
#' # Plot R-squared values without phylogeny
#' plot(result, what = "rsquared", plot_phylogeny = FALSE)
#'
#' # Plot a PCA of coefficients with superimposed the phylogeny
#' # (a so called "phylomorphospace" plot)
#' plot(result, what = "coefficients")
#'
#' # Plot coefficients without phylogeny
#' plot(result, what = "coefficients", plot_phylogeny = FALSE)
#'
#' # Single trait example
#' Y_single = rnorm(20)
#' result_single = phylo_w_regr_fit(tree, Y_single, X)
#' plot(result_single, what = "coefficients")
#'
#' @export
phylo_w_regr_fit = function(tree, Y, X) {
  # Compute the phylogenetic correlation matrix (weights) using Brownian motion model
  wgt = nlme::corMatrix(nlme::Initialize(ape::corBrownian(phy=tree, form= ~tree$tip.label), as.data.frame(Y)))


  # For each tip in the tree, fit a weighted linear model of Y ~ X
  PWR = lapply(tree$tip.label, function(x) {
    # Fit weighted linear regression for tip x
    model_w = stats::lm(Y ~ X, weights=wgt[,x], y=TRUE)
    # Calculate R-squared for the fitted model
    if (is.null(dim(Y)) || ncol(Y)==1) {
      rsquared = var(model_w$fitted.values) / var(model_w$y)
    } else {
      rsquared = sum(apply(model_w$fitted.values, 2, var)) /
        sum(apply(model_w$y, 2, var))
    }
    # Extract coefficients and intercept
    # Handle both vector and matrix cases for coefficients
    if (is.matrix(model_w$coefficients)) {
      # Matrix case: multiple dependent variables
      coefficients = model_w$coefficients[-1,] # Exclude intercept
      intercept = model_w$coefficients[1,]     # Intercept only
    } else {
      # Vector case: single dependent variable
      coefficients = model_w$coefficients[-1]  # Exclude intercept
      intercept = model_w$coefficients[1]      # Intercept only
    }

    results = list(
      rsquared = rsquared,
      coefficients = coefficients,
      intercept = intercept
    )
    return(results)
  })

  # Assign tip labels as names for the results
  names(PWR) = tree$tip.label

  # Collect R-squared, intercepts, and coefficients for all tips
    rsquared = sapply(PWR, function(x) x$rsquared)
    intercepts = t(sapply(PWR, function(x) x$intercept))
    coeff = t(sapply(PWR, function(x) x$coefficients))
    # The part below handles the univariate case
    # where intercepts and coefficients are vectors
    if(nrow(intercepts)==1){
      intercepts=as.vector(intercepts)
      coeff=as.vector(coeff)
      names(intercepts)=names(rsquared)
      names(coeff)=names(rsquared)
    }

  # Return a list containing all results and the weights matrix
  results=list(
    tree = tree,
    rsquared = rsquared,
    intercepts = intercepts,
    coefficients = coeff,
    weights = wgt
  )
  class(results)=c("phylo_wregr_fit","list")
 return(results)
}

#' Plot method for phylo_wregr_fit objects
#'
#' Visualizes the results of phylogenetic weighted regression analysis. The function can plot
#' either R-squared values or regression coefficients, with or without phylogenetic context.
#' For multivariate coefficients, Principal Component Analysis is performed and the first two
#' components are plotted. When phylogenetic context is included, the plots are overlaid on
#' the phylogenetic tree structure.
#'
#' @param x An object of class `phylo_wregr_fit`, typically the result of \code{\link{phylo_w_regr_fit}}
#' @param what Character string specifying what to plot. Either "rsquared" for R-squared values
#'   or "coefficients" for regression coefficients. Default is "rsquared".
#' @param plot_phylogeny Logical indicating whether to include phylogenetic information in plots.
#'   When TRUE, uses phylogenetic plotting functions from the phytools package. Default is TRUE.
#' @param ... Additional arguments passed to the underlying plotting functions (e.g., graphical parameters)
#'
#' @details
#' The function uses different visualization approaches depending on the data type and user preferences:
#'
#' \strong{For R-squared values:}
#' \itemize{
#'   \item With phylogeny: Uses \code{phytools::plotTree.wBars} to show values as bars along tree tips
#'   \item Without phylogeny: Uses \code{stats::density} to show the distribution of R-squared values
#' }
#'
#' \strong{For coefficients:}
#' \itemize{
#'   \item Multivariate case with phylogeny: Performs PCA using \code{stats::prcomp} and plots the first
#'     two components in phylomorphospace using \code{phytools::phylomorphospace}
#'   \item Multivariate case without phylogeny: Plots PC1 vs PC2 scores using base \code{plot}
#'   \item Univariate case with phylogeny: Uses \code{phytools::plotTree.wBars}
#'   \item Univariate case without phylogeny: Uses \code{stats::density}
#' }
#'
#' @return Invisibly returns the input object x
#'
#' @note The phytools package is required when \code{plot_phylogeny = TRUE}. The function will
#' throw an error if phytools is not installed and phylogenetic plots are requested.
#'
#' @seealso \code{\link{phylo_w_regr_fit}} for the main function that generates objects of this class
#'
#' @method plot phylo_wregr_fit
#' @importFrom stats prcomp density
#'
#' @examples
#' # Same example as the main function but with more complex plotting
#' # using the attributes passed to functions in phytools
#' library(phytools)
#' library(viridisLite)
#'
#' set.seed(123)
#'
#' tree = pbtree(n=20)
#' # Generate a tre with 20 tips
#'
#' Y = matrix(rnorm(100), nrow=20)
#' X = rnorm(20)
#' result = phylo_w_regr_fit(tree, Y, X)
#'
#' color_inferno=inferno(n=100, direction = -1)
#' # Create a color palette for the R-squared values using 100 values
#' # with lighter colours representing smaller values
#' # see also ?inferno
#' color_result_rsquared=sapply(result$rsquared,
#'                             function(x) color_inferno[round(x*100)])
#' # Map the values of rsquared to the color palette
#' plot(result, what = "rsquared", plot_phylogeny = TRUE,
#'      type = "phylogram", col = color_result_rsquared)
#' # Plot the phylogeny with the bars and the colours representing the R-squared values
#' # Now add a legend of colours and rsquared value covering all range in color_turbo
#'
#'
#'
#'
#' @export
plot.phylo_wregr_fit = function(x, what = c("rsquared", "coefficients"), plot_phylogeny = TRUE, ...) {
  what = match.arg(what)

  # Check for required package if phylogeny plots are requested
  if (plot_phylogeny && !requireNamespace("phytools", quietly = TRUE)) {
    stop("Package 'phytools' is required for plots involving the phylogeny. Please install it.")
  }

  if (what == "coefficients") {
    coefs = x$coefficients
    if (is.matrix(coefs) && ncol(coefs) > 1) {
      # PCA on covariance matrix of coefficients
      pca = prcomp(coefs, center = TRUE, scale. = FALSE)
      scores = pca$x
      if (plot_phylogeny) {
        # Phylomorphospace plot
        phytools::phylomorphospace(x$tree, scores[, 1:2], label = "horizontal", ...)
      } else {
        # Scatterplot of first two PC scores
        plot(scores[, 1], scores[, 2], asp=1, xlab = "PC1", ylab = "PC2", main = "PCA of Coefficients", ...)
      }
    } else {
      # Vector case
      vals = as.vector(coefs)
      if (plot_phylogeny) {
        vals = setNames(as.vector(vals), x$tree$tip.label)
        stopifnot(all(names(vals) == x$tree$tip.label))
        phytools::plotTree.wBars(x$tree, vals, ...)
      } else {
        plot(density(vals), main = "Density of Coefficients", xlab = "Coefficient", ...)
      }
    }
  } else if (what == "rsquared") {
    vals = as.vector(x$rsquared)
    if (plot_phylogeny) {
      vals = setNames(as.vector(vals), x$tree$tip.label)
      stopifnot(all(names(vals) == x$tree$tip.label))
      phytools::plotTree.wBars(x$tree, vals, ...)
    } else {
      plot(density(vals), main = "Density of R-squared", xlab = "R-squared", ...)
    }
  }

  invisible(x)
}
