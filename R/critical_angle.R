#' Compute the critical angle for the test of the angle between two
#' multivariate vectors
#'
#' This function computes the angle below which two vectors would be
#' "significantly similar" using Li 2011 test.
#'
#' This function considers the formulas in Li 2011 which have been
#' used to perform a test of the angle between multivariate vectors.
#' This is the test implemented in MorphoJ (Klingenberg 2011),
#' in the R package Morpho (Schlager 2016), and in the function
#' TestOfAngle of this package.
#' The test produces a p value for the angle between two vectors
#' (with significant values being "significantly similar").
#'
#' This function reverts the logic of the formula/test and, given a
#' number of dimensions (e.g., morphometric variables) and a level
#' of significance (by default 0.05),
#' it computes the angle below which the test would be significant.
#'
#'
#'
#' @section Notice:
#' In case of very many dimensions, there are numerical problems in
#' performing the test.
#' If prec="normal" the function uses internally the package Morpho
#' (which should be installed) and these problems will start occurring
#' above about 340 variables.
#' If prec="mpfr" the function uses internally a custom function which
#' requires the package Rmpfr (which should be installed).
#' This allows the computation of extremely large or small numbers, it
#' is currently set to a 120 bit precision and allows the computation
#' using up to about 1200 variables (over this number, there will be
#' problems even using the option "mpfr")
#'
#' @param dimensions number of dimensions to use
#' @param alpha significance level (by default 0.05)
#' @param prec whether one has to use the internal R functions
#'   ("normal", default), or mpfr precision ("mpfr")
#'
#' @return The function outputs the critical angle (in degrees)
#'
#' @references Klingenberg CP. 2011. MorphoJ: an integrated software
#' package for geometric morphometrics. Mol Ecol Resour 11:353-357.
#' @references Li S. 2011. Concise formulas for the area and volume of
#' a hyperspherical cap. Asian Journal of Mathematics and Statistics
#' 4:66-70.
#' @references Schlager S. 2016. Morpho: Calculations and
#' Visualisations Related to Geometric Morphometrics.
#'
#' @seealso \code{\link{TestOfAngle}}
#'
#' @examples
#'
#' critical_angle(50)
#' # This is the angle below which two vectors (e.g., PCA axes)
#' # of 50 variables would be would be "significantly similar"
#'
#' @export
critical_angle=function(dimensions, alpha=0.05, prec="normal") {

  # Validate inputs
  if (!is.numeric(dimensions) || dimensions <= 0 || dimensions %% 1 != 0) {
    stop("dimensions must be a positive integer")
  }
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be between 0 and 1")
  }
  if (!prec %in% c("normal", "mpfr")) {
    stop("prec must be either 'normal' or 'mpfr'")
  }


  n=dimensions
  if(prec=="mpfr") {
  SP=areaSphere_mpfr(n)
  } else {
    SP=Morpho::areaSphere(n)
  }
  # Area of the hypersphere

  Apart_target=SP*alpha
  ib=Apart_target/(0.5*SP)
  QB=qbeta(as.numeric(ib), (n-1)/2,0.5)
  ASIN=asin(sqrt(QB))
return(rad2deg(ASIN))
}

areaSphere_mpfr <- function(n,r=1) {
  nom <- 2*pi^(n/2)*r^(n-1)
  denom <- gamma(Rmpfr::mpfr(n/2, precBits=120))
return(nom/denom)
}


