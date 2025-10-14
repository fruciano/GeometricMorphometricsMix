#' Permutation test of difference in disparity/morphospace occupation
#'
#' Performs a permutation test of difference in disparity between two
#' groups.
#'
#' The function employs commonly used test statistics
#' to quantify disparity/morphospace occupation/variation in each group.
#' The two statistics currently implemented are multivariate variance
#' (also known as sum of variances, trace of the covariance matrix,
#' Procrustes variance),
#' and mean pairwise Euclidean distances.
#' These two metrics have a long history in the quantification of
#' disparity both in geometric morphometrics
#' (e.g., Zelditch et al. 2004; Fruciano et al., 2014, 2016) and more
#' in general in evolution (e.g., Foote, 1996; Willis 2001)
#' The observed statistics are then compared to their empirical
#' distributions obtained through permutations,
#' to obtain a p-value.
#'
#' @section Notice:
#' The values of the test statistics in the output are the observed
#' in the sample.
#' If they are of interest, and the two groups have different sample
#' size, consider computing their rarefied versions (for instance with
#' the function \code{\link{disparity_resample}})
#' for reporting in papers and the like.
#'
#' @param X1,X2 Matrices or data frames containing data for each group
#' (observations in rows, variables in columns).
#' @param perm number of permutations
#' @seealso \code{\link{disparity_resample}}, \code{\link{BTailTest}}
#'
#' @return The function outputs a dataframe containing:
#'   the observed values of the tests statistics for each group,
#'   their absolute differences, and
#'   The p values obtained through the permutational procedure
#'
#' @references Foote M. 1997. The evolution of morphological diversity.
#' Annual Review of Ecology and Systematics 28:129.
#' @references Wills MA. 2001. Morphological disparity: a primer. In.
#' Fossils, phylogeny, and form: Springer. p. 55-144.
#' @references Zelditch ML, Swiderski DL, Sheets HD. 2004. Geometric
#' morphometrics for biologists: a primer: Academic Press.
#' @references Fruciano C, Franchini P, Raffini F, Fan S, Meyer A.
#' 2016. Are sympatrically speciating Midas cichlid fish special?
#' Patterns of morphological and genetic variation in the closely
#' related species Archocentrus centrarchus. Ecology and Evolution
#' 6:4102-4114.
#' @references Fruciano C, Pappalardo AM, Tigano C, Ferrito V. 2014.
#' Phylogeographical relationships of Sicilian brown trout and the
#' effects of genetic introgression on morphospace occupation.
#' Biological Journal of the Linnean Society 112:387-398.
#'
#' @examples
#' library(MASS)
#' set.seed(123)
#'
#' X1=mvrnorm(20, mu=rep(0, 40), Sigma=diag(40))
#' X2=mvrnorm(100, mu=rep(5, 40), Sigma=diag(40))
#' # create two groups of random observations
#' # with different means and sample sizes,
#' # but the same covariance matrix
#'
#' # We expect that the two groups will have the same
#' # variance (disparity/morphospace occupation)
#' # and therefore the test will be non-significant
#'
#' disparity_test(X1, X2, perm=999)
#' # This is, indeed, the case
#'
#' @export
disparity_test=function(X1, X2, perm=999) {

  observed_muvar=c(muvar(X1), muvar(X2))
  observed_mpED=c(
    meanpairwiseEuclideanD(X1), meanpairwiseEuclideanD(X2)
  )
  diff_muvar=abs(observed_muvar[1]-observed_muvar[2])
  diff_mpED=abs(observed_mpED[1]-observed_mpED[2])
  # Observed statistics and their differences

  n_grp1=nrow(X1)
  n_tot=n_grp1+nrow(X2)
  n_tot_idx=seq_len(n_tot)

  X1X2=rbind(X1, X2)
  PermutedDatasets=lapply(seq_len(perm), function(x) {
      X1X2perm=X1X2[sample(n_tot_idx, n_tot, replace=FALSE),]
        list(X1perm=X1X2perm[seq_len(n_grp1),],
             X2perm=X1X2perm[(n_grp1+1):n_tot,])
            }
          )
  perm_diff_muvar=unlist(lapply(PermutedDatasets, function(x)
    abs(muvar(x[[1]])-muvar(x[[2]])))
  )
  perm_diff_mpED=unlist(lapply(PermutedDatasets, function(x)
    abs(meanpairwiseEuclideanD(x[[1]])-
          meanpairwiseEuclideanD(x[[2]])))
  )
  # Create permuted datasets and compute
  # the absolute difference in disparity estimates
  # between permuted samples

  p_value_muvar=(length(which(perm_diff_muvar >= diff_muvar))+1)/(perm + 1)
  p_value_mpED=(length(which(perm_diff_mpED >= diff_mpED))+1)/(perm + 1)
  # Compute p_values

  Results=data.frame( cbind(rbind(observed_muvar, observed_mpED),
                      difference=c(diff_muvar, diff_mpED),
                      p_value=c(p_value_muvar, p_value_mpED)
                        )
                    )
  colnames(Results)[1:2]=c("Observed_grp1", "Observed_grp2")
  rownames(Results)=c(
    "Multivariate variance", "Mean pairwise Euclidean distance"
  )
  # Create a table of results

return(Results)
}




