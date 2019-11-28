#' Rarefied estimates of disparity/morphospace occupation
#'
#' Computes rarefied estimates of multivariate variance
#' (i.e., the trace of the covariance matrix) and
#' the mean pairwise Euclidean distance
#'
#' The advantage of computing rarefied estimates of a parameter
#' is to account for differences in sample sizes.
#' This can be useful to compare estimates among groups.
#'
#' This function is similar to the approaches implmented in the
#' Matlab package MDA (Navarro 2003).
#'
#' Examples of these approaches abound in the literature
#' (e.g., Drake & Klingenberg 2010; Fruciano et al. 2014, 2016)
#'
#' Optionally, the function also computes "Proper Variance" (Claramunt 2010)
#' using linear shrinkage estimates of the covariance matrix.
#'
#' @section Notice:
#' the function does NOT perform GPA on each rarefied sample
#' this may or may not make a difference in estimates.
#' In many cases, it will probably not make much difference
#'
#' @section Notice:
#' If claramunt_pv=TRUE the package nlshrink should be installed.
#' As this computation is based on linear shrinkage estimates of the covariance matrix
#' it can be computational intensive.
#' Please, also insure to check (and cite) the original publication (Claramunt 2010)
#' to understand better this estimator.
#'
#' @param Data Matrix or data frame containing data
#' (observations in rows, variables in columns).
#' @param rep number of resamplings to obtain the rarefied estimate
#' @param samplesize sample size to which the rarefaction procedure is carried out
#' @param claramunt_pv if TRUE, compute an estimate of "Proper variance" (Claramunt 2010)
#' @seealso \code{\link{rarefied_convex_hull}}, \code{\link{BTailTest}}
#'
#' @return The function outputs a list with the following elements:
#'  \describe{
#'   \item{RarefiedSamplesEstimates}{Estimates for each of the rarefied samples}
#'   \item{Descriptives}{Mean, standard deviation, median and 2.5\% (Min) and 97.5\% (Max) across the rarefied samples}
#' }
#'
#' @references Navarro N. 2003. MDA: a MATLAB-based program for morphospace-disparity analysis. Computers & Geosciences 29:655-664.
#' @references Drake AG, Klingenberg CP. 2010. Large-scale diversification of skull shape in domestic dogs: disparity and modularity. American Naturalist 175:289-301.
#' @references Claramunt S. 2010. Discovering exceptional diversifications at continental scales: the case of the endemic families of Neotropical Suboscine passerines. Evolution 64:2004-2019.
#' @references Fruciano C, Franchini P, Raffini F, Fan S, Meyer A. 2016. Are sympatrically speciating Midas cichlid fish special? Patterns of morphological and genetic variation in the closely related species Archocentrus centrarchus. Ecology and Evolution 6:4102-4114.
#' @references Fruciano C, Pappalardo AM, Tigano C, Ferrito V. 2014. Phylogeographical relationships of Sicilian brown trout and the effects of genetic introgression on morphospace occupation. Biological Journal of the Linnean Society 112:387-398.
#'
#' @import stats
#' @export
rarefied_disparity=function(Data,rep=1000,samplesize, claramunt_pv=FALSE) {

    RarefiedSamples=lapply(seq(rep),function(x) Data[sample(1:nrow(Data),samplesize,replace=TRUE),])

    if (claramunt_pv==TRUE) {
      DisparityEst=lapply(RarefiedSamples, function(x)
          c(muvar(x),meanpairwiseEuclideanD(x), claramunt_proper_variance(x)))
      DisparityEst=do.call(rbind,DisparityEst)
      colnames(DisparityEst)=c("Multivariate variance",
                               "Mean pairwise Euclidean distance",
                               "Claramunt proper variance")

      Descriptives=apply(DisparityEst, 2, function (x)
        c(Mean=mean(x), SD=sd(x), Median=quantile(x,0.5), Min=quantile(x,0.025), Max=quantile(x,0.975) )
      )
      Results=list(RarefiedSamplesEstimates=DisparityEst,Descriptives=Descriptives)
    } else {

    DisparityEst=lapply(RarefiedSamples, function(x) c(muvar(x),meanpairwiseEuclideanD(x)))
    DisparityEst=do.call(rbind,DisparityEst)
    colnames(DisparityEst)=c("Multivariate variance","Mean pairwise Euclidean distance")

    Descriptives=apply(DisparityEst, 2, function (x)
    c(Mean=mean(x), SD=sd(x), Median=quantile(x,0.5), Min=quantile(x,0.025), Max=quantile(x,0.975) )
    )
    Results=list(RarefiedSamplesEstimates=DisparityEst,Descriptives=Descriptives)
    }
return(Results)
}
