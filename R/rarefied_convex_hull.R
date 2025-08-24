#' Rarefied estimates of n-dimensional convex hull volume
#'
#' @description
#' 
#' `rarefied_convex_hull()` is deprecated and will be removed in a future version.
#'
#' Computes rarefied estimates of convex hull volume
#' (i.e., the n-dimensional volume of the convex hull containing
#' a set of observations).
#'
#' The advantage of computing rarefied estimates of a parameter
#' is to account for differences in sample sizes.
#' This can be useful to compare estimates among groups.
#'
#' n-dimensional convex hulls are not frequently used in geometric morphometrics
#' to quantify morphospace occupation
#' (but see Drake & Klingenberg 2010; Fruciano et al. 2012, 2014).
#' However, they can be useful to quantify the spread of the points
#' (observations, individuals) in multivariate (tangent shape) space.
#' Because of their own nature, they are obviously sensitive to outliers.
#' For this reason, rarefaction - as implemented here - is particularly useful
#' to compare patterns across groups.
#'
#' @section Notice:
#' The function uses internally the package geometry, which should be installed
#' prior to using this function.
#' Notice also that computation of multidimensional convex hulls is computationally expensive
#' and that requires more observations (individuals) than variables.
#' For this reason, errors might be produced when using too small sample sizes and/or too many variables.
#' A solution for this is increasing sample sizes or performing dimensionality reduction
#' by using this approach only on a meaningful subset of principal components.
#' The function does NOT perform GPA on each rarefied sample
#' this may or may not make a difference in estimates.
#' In many cases, it will probably not make much difference.
#'
#' @section Citation:
#' If you use this function, in addition to this package,
#' please cite Fruciano et. al (2014), which is the first use of
#' rarefied n-dimensional convex hull in morphometrics.
#'
#' @param Data Matrix or data frame containing data
#' (observations in rows, variables in columns).
#' @param rep number of resamplings to obtain the rarefied estimate
#' @param samplesize sample size to which the rarefaction procedure is carried out
#' @seealso \code{\link{rarefied_disparity}}, \code{\link{BTailTest}}
#'
#' @return The function outputs a list with the following elements:
#'  \describe{
#'   \item{RarefiedSamplesEstimates}{Estimates of both multivariate variance and
#'         mean pairwise Euclidean distance for each of the rarefied samples}
#'   \item{Descriptives}{Mean, standard deviation, median and 2.5\% (Min) and 97.5\% (Max) across the rarefied samples}
#' }
#'
#' @references Drake AG, Klingenberg CP. 2010. Large-scale diversification of skull shape in domestic dogs: disparity and modularity. American Naturalist 175:289-301.
#' @references Fruciano C, Tigano C, Ferrito V. 2012. Body shape variation and colour change during growth in a protogynous fish. Environmental Biology of Fishes 94:615-622.
#' @references Fruciano C, Pappalardo AM, Tigano C, Ferrito V. 2014. Phylogeographical relationships of Sicilian brown trout and the effects of genetic introgression on morphospace occupation. Biological Journal of the Linnean Society 112:387-398.
#'
#' @import stats
#' @export
rarefied_convex_hull=function(Data,rep=100,samplesize) {
    .Deprecated(msg = "rarefied_convex_hull() is deprecated and will be removed in a future version. Please, use the more general disparity_resample() instead.")
    
    RarefiedSamples=lapply(seq(rep),function(x) Data[sample(1:nrow(Data),samplesize,replace=TRUE),])

    DisparityEst=lapply(RarefiedSamples, function(x) geometry::convhulln(x,options="FA")$vol)
    DisparityEst=do.call(rbind,DisparityEst)
    colnames(DisparityEst)="Convex hull volume"

    Descriptives=apply(DisparityEst, 2, function (x)
    c(Mean=mean(x), SD=sd(x), Median=quantile(x,0.5), Min=quantile(x,0.025), Max=quantile(x,0.975) )
    )
    Results=list(RarefiedSamplesEstimates=DisparityEst,Descriptives=Descriptives)
return(Results)
}


