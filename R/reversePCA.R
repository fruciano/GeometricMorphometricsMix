#' 'Reverse' PCA
#'
#' Obtain data in the original variables starting from PCA scores
#'
#' Given a set of PCA scores, a set of eigenvectors and the mean of the data on which
#' the PCA was originally computed, this function returns a set of observations in the
#' original data space.
#' This simple function can have many applications (not only in morphometrics) such as
#'  \itemize{
#'   \item Producing predicted shapes for extreme, unobserved, values of PC scores
#'   (e.g., the shape, or other set of variables, corresponding to 3 times the most
#'   extreme positive PC1 score), or their combinations (e.g., a value of 0.4 on PC1
#'   and 0.3 on PC2)
#'   \item Interpret results of analyses carried out on PC scores by converting these
#'   back to shapes
#' }
#'
#' @param Scores matrix n x s of n observation for s scores. Can be a single column
#' (for instance, if one is interested in PC1 only)
#' @param Eigenvectors matrix p x s of eigenvectors (as column eigenvectors);
#' notice that s (number of column eigenvectors) must be the same as the number of columns
#' of PC scores in Scores
#' @param Mean vector of length p with the multivariate mean computed on the original dataset
#' (the dataset on which the PCA was carried out)
#'
#' @return The function outputs a matrix n x p of the scores 'back-transformed'
#' into the original variables
#'
#' @examples
#'
#' library(MASS)
#' set.seed(123)
#' A=mvrnorm(10,mu=rep(1,2),Sigma=diag(2))
#' # Create a small dataset of 10 observations and 2 variables
#'
#' Mean=colMeans(A)
#' PCA=prcomp(A)
#' # Compute mean and perform PCA
#'
#' B=reversePCA(Scores=PCA$x, Eigenvectors=PCA$rotation, Mean=Mean)
#' # 'Revert' the PCA (in this case using all scores and all axes
#' # to get the same results as the starting data)
#'
#' round(A,10)==round(B,10)
#' # Check that all the results are the same
#' # (rounding to the 10th decimal because of numerical precision)
#'
#'
#' @import stats
#' @export
reversePCA = function(Scores, Eigenvectors, Mean) {
    ZEt = Scores %*% t(Eigenvectors)
    # Multiply the scores by the transpose of the eigenvectors
    rows = nrow(Scores)
    X = ZEt + rep.row(Mean, rows)
    # Add to the mean
    return(X)
}
