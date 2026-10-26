#' Permutation test of difference in disparity/morphospace occupation
#'
#' Performs a permutation test of difference in disparity between two
#' groups using multivariate variance and mean pairwise Euclidean distance.
#'
#' @details
#' This compatibility interface calls \code{\link{disparity_compare}} with
#' both disparity measures, absolute differences, and pooled permutations.
#' The arguments and output retain the interface of earlier package versions.
#' Both statistics share the same permutations. P-values use the upper
#' empirical tail, including ties, with the plus-one correction.
#' Pooled permutation requires exchangeability under the null; equality of
#' disparity alone does not establish this assumption 
#' (see also \code{\link{disparity_compare}} for more details).
#'
#' @section Notice:
#' Reported disparity estimates are observed, rather than rarefied, values.
#' Consider \code{\link{disparity_resample}} when sample-size standardization
#' is needed for reporting disparity estimates.
#'
#' @param X1,X2 Numeric matrices or data frames containing data for each
#'   group (observations in rows, variables in columns), with matching
#'   variables in the same order. Each group must have at least two rows.
#' @param perm Positive integer number of permutations. Default is 999.
#'
#' @return A data frame with columns \code{Observed_grp1},
#'   \code{Observed_grp2}, \code{difference} (absolute difference), and
#'   \code{p_value}. Rows are \code{"Multivariate variance"} and
#'   \code{"Mean pairwise Euclidean distance"}, in that order.
#'
#' @seealso \code{\link{disparity_compare}},
#'   \code{\link{disparity_resample}}, \code{\link{BTailTest}}
#' @references Foote M. 1997. The evolution of morphological diversity.
#'   Annual Review of Ecology and Systematics 28:129-152.
#' @references Wills MA. 2001. Morphological disparity: a primer. In:
#'   Fossils, phylogeny, and form. Springer. p. 55-144.
#' @references Zelditch ML, Swiderski DL, Sheets HD. 2004. Geometric
#'   morphometrics for biologists: a primer. Academic Press.
#' @references Fruciano C, Franchini P, Raffini F, Fan S, Meyer A. 2016.
#'   Are sympatrically speciating Midas cichlid fish special? Patterns of
#'   morphological and genetic variation in the closely related species
#'   Archocentrus centrarchus. Ecology and Evolution 6:4102-4114.
#' @references Fruciano C, Pappalardo AM, Tigano C, Ferrito V. 2014.
#'   Phylogeographical relationships of Sicilian brown trout and the effects
#'   of genetic introgression on morphospace occupation.
#'   Biological Journal of the Linnean Society 112:387-398.
#'
#' @examples
#' set.seed(123)
#' X1=matrix(rnorm(80), nrow=20)
#' X2=matrix(rnorm(140), nrow=35)
#' disparity_test(X1, X2, perm=99)
#'
#' @export
disparity_test=function(X1, X2, perm=999) {
  if (!(is.matrix(X1) || is.data.frame(X1)) ||
      !(is.matrix(X2) || is.data.frame(X2)) || ncol(X1) != ncol(X2)) {
    stop("'X1' and 'X2' must be matrices or data frames with matching variables.")
  }
  comparison=disparity_compare(
    Data=rbind(X1, X2),
    group=factor(rep(c("grp1", "grp2"), c(nrow(X1), nrow(X2))),
                 levels=c("grp1", "grp2")),
    statistic=c("multivariate_variance", "mean_pairwise_euclidean_distance"),
    comparison="absolute_difference", method="permutation", n_resamples=perm
  )
  Results=comparison$results[,c("Observed_grp1", "Observed_grp2",
                                "contrast", "p_value")]
  colnames(Results)[3]="difference"
  return(Results)
}
