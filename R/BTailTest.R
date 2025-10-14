#' BTailTest for difference in disparity/morphospace occupation
#'
#' Performs the BTailTest, in the same spirit as implemented in the
#' Matlab package MDA (Navarro 2003) and used in various empirical
#' papers (e.g., Fruciano et al. 2014, 2016).
#'
#' This is a test of the difference in disparity between two groups:
#' a reference and a test group.
#' The function proceeds by computing a bootstrapped distribution of
#' the test statistics (multivariate variance and mean pairwise
#' Euclidean distances in this implementation) in the reference sample
#' and then comparing the statistics observed in the test sample to
#' this distribution to obtain p-values.
#'
#'
#' @section Citation:
#' If you use this function, in addition to this package,
#' please cite Navarro (2003)
#'
#' @param Reference Matrix or data frame containing data for the
#'   reference group (observations in rows, variables in columns).
#' @param Test Matrix or data frame containing data for the test
#'   group (observations in rows, variables in columns).
#' @param boot number of bootstrap replicates
#' @seealso \code{\link{disparity_resample}},
#'   \code{\link{disparity_test}}
#'
#' @return The function outputs a list with the following elements:
#'  \describe{
#'   \item{BootstrappedSamplesEstimates}{Estimates of both
#'         multivariate variance and mean pairwise Euclidean distance
#'         for each of the bootstrapped samples}
#'   \item{pvalues}{p values obtained for the test}
#' }
#'
#' @references Navarro N. 2003. MDA: a MATLAB-based program for
#'   morphospace-disparity analysis. Computers & Geosciences
#'   29:655-664.
#' @references Fruciano C, Franchini P, Raffini F, Fan S, Meyer A.
#'   2016. Are sympatrically speciating Midas cichlid fish special?
#'   Patterns of morphological and genetic variation in the closely
#'   related species Archocentrus centrarchus. Ecology and Evolution
#'   6:4102-4114.
#' @references Fruciano C, Pappalardo AM, Tigano C, Ferrito V. 2014.
#'   Phylogeographical relationships of Sicilian brown trout and the
#'   effects of genetic introgression on morphospace occupation.
#'   Biological Journal of the Linnean Society 112:387-398.
#' @export
BTailTest=function(Reference,Test,boot=1000) {

	Reference=scale(Reference, center=TRUE, scale=FALSE)
	Test=scale(Test, center=TRUE, scale=FALSE)

ObservedTestStats=c(muvar(Test),meanpairwiseEuclideanD(Test))

BootstrappedSamples=lapply(seq(boot), function(x)
  Reference[sample(seq_len(nrow(Reference)), nrow(Reference),
                   replace=TRUE), ])

DisparityEst=lapply(BootstrappedSamples, function(x)
  c(muvar(x), meanpairwiseEuclideanD(x)))
DisparityEst=do.call(rbind,DisparityEst)
colnames(DisparityEst)=c("Multivariate variance",
                         "Mean pairwise Euclidean distance")

pvalues=lapply(seq_along(ObservedTestStats), function(x) c(
  length(which(DisparityEst[,x]>ObservedTestStats[x]))/boot,
  length(which(DisparityEst[,x]<ObservedTestStats[x]))/boot,
  2*min(length(which(DisparityEst[,x]>ObservedTestStats[x]))/boot,
        length(which(DisparityEst[,x]<ObservedTestStats[x]))/boot)
))

pvalues=do.call(rbind,pvalues)
colnames(pvalues)=c("Upper tail","Lower tail","Two-tail")
rownames(pvalues)=c("Multivariate variance",
                    "Mean pairwise Euclidean distance")

Results=list(BootstrappedSamplesEstimates=DisparityEst,pvalues=pvalues)
return(Results)
}

