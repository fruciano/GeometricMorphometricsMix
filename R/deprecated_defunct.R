#' Defunct functions in GeometricMorphometricsMix
#'
#' @description
#' These functions are defunct and have been removed from the package.
#' Please use the modern alternatives as indicated below.
#'
#' @name GeometricMorphometricsMix-defunct
#' @keywords internal
NULL

#' @rdname GeometricMorphometricsMix-defunct
#' @section \code{rarefied_convex_hull}:
#' For \code{rarefied_convex_hull}, use \code{\link{disparity_resample}} with 
#' \code{statistic='convex_hull_volume'} instead.
#' @export
rarefied_convex_hull = function(Data, rep = 100, samplesize) {
    .Defunct(msg = "rarefied_convex_hull() is defunct. Please use disparity_resample() with statistic='convex_hull_volume' instead.")
}

#' @rdname GeometricMorphometricsMix-defunct  
#' @section \code{rarefied_disparity}:
#' For \code{rarefied_disparity}, use \code{\link{disparity_resample}} with 
#' appropriate statistic parameter instead.
#' @export
rarefied_disparity = function(Data, rep = 1000, samplesize, claramunt_pv = FALSE) {
    .Defunct(msg = "rarefied_disparity() is defunct. Please use disparity_resample() with appropriate statistic parameter instead.")
}
