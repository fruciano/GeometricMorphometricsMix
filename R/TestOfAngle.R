#' Perform a test of the angle between two multivariate vectors
#'
#' This function performs a test for the angle between two
#' multivariate vector, optionally allowing for "flipping" one of
#' the axes
#'
#' This function is useful to perform a test of the angle between
#' two multivariate vectors (e.g., principal component eigenvectors
#' computed from two covariance matrices, such as covariance matrices
#' of two different species).
#' The function uses internally angleTest from the package Morpho by
#' Stefan Schlager.
#' That function, in turn, uses the formulas in Li 2011 to compute
#' significance, rather than using permutations.
#' This is the same approach also implemented in MorphoJ
#' (Klingenberg 2011).
#' There is no special advantage in using this function relative to
#' the original in the package Morpho, except that this function
#' outputs angles both in radians and degrees and that it optionally
#' allow to "flip" one of the two axes. This is useful in the cases
#' where the direction of one axis is arbitrary, such as in PCA
#' (where positive and negative values are interchangable and
#' recomputing PCA might result in positive scores for the
#' observations which were negative (and vice versa).
#' In this situation, angles larger than 90 degrees are not
#' meaningful and one way of dealing with this is, by choosing the
#' option flip=TRUE (default), to change the sign of one of the two
#' vectors ("flip").
#'
#' @section Notice:
#' The p value is for the probability that the angle between two
#' random vectors is smaller or equal to the one computed from the
#' vectors provided (x, y). This means that significant values
#' indicate that the two provided vectors are "significantly
#' similar", whereas non-significant values means that the two
#' vectors are substantially different
#'
#' @param x,y numerical vectors
#' @param flip logical stating whether (TRUE, default) axes should
#'   be "flipped" in case the angle is larger than 90 degrees
#'   (see Details)
#'
#' @return The function outputs a list with
#' \item{angle}{angle between vectors in radians}
#' \item{angle_deg}{angle between vectors in degrees}
#' \item{p.value}{p value}
#' \item{critical_angle}{angle below which two vectors of the same
#'   number of dimensions as the ones being tested would be
#'   significant}
#'
#' @references Klingenberg CP. 2011. MorphoJ: an integrated software
#'   package for geometric morphometrics. Mol Ecol Resour
#'   11:353-357.
#' @references Li S. 2011. Concise formulas for the area and volume
#'   of a hyperspherical cap. Asian Journal of Mathematics and
#'   Statistics 4:66-70.
#' @references Schlager S. 2016. Morpho: Calculations and
#'   Visualisations Related to Geometric Morphometrics.
#'
#' @seealso \code{\link{critical_angle}}
#'
#' @examples
#' library(MASS)
#' library(clusterGeneration)
#' set.seed(123)
#' Data=mvrnorm(500,mu=rep(1,50),Sigma=genPositiveDefMat(dim=50)$Sigma)
#' A=Data[1:250,]
#' B=Data[251:500,]
#' # Create two groups of observations (e.g., specimens)
#' # from the same multivariate normal distribution
#' # with the same starting covariance matrix
#'
#' PCA_A=prcomp(A)
#' PCA_B=prcomp(B)
#' # Perform separate PCAs for each of the two groups of observations
#'
#' TestOfAngle(PCA_A$rotation[,1],PCA_B$rotation[,1], flip=TRUE)
#' # Test for the angle between the two first eigenvectors
#'
#' @export
TestOfAngle=function(x,y, flip=TRUE) {
AngleTest=Morpho::angleTest(x,y)
if (all(c(flip, AngleTest$angle > pi/2))) {
AngleTest=Morpho::angleTest((-1*x),y)
}
AngleTest$angle_deg=rad2deg(AngleTest$angle)
AngleTest$critical_angle=critical_angle(length(x))
return(AngleTest)
}



