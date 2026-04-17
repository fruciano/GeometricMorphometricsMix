#' Parallel implementation of Adams' Kmult (defunct)
#'
#' @description
#' **Defunct**: This function has been removed from the package.
#' Please use \code{geomorph::physignal()} instead.
#'
#' @param data Ignored (function is defunct).
#' @param trees Ignored (function is defunct).
#' @param burninpercent Ignored (function is defunct).
#' @param iter Ignored (function is defunct).
#' @param verbose Ignored (function is defunct).
#'
#' @return Does not return; always signals a defunct error.
#'
#' @references Adams DC. 2014. A Generalized K Statistic for Estimating
#' Phylogenetic Signal from Shape and Other High-Dimensional
#' Multivariate Data. Systematic Biology 63:685-697.
#' @references Fruciano C, Celik MA, Butler K, Dooley T, Weisbecker V,
#' Phillips MJ. 2017. Sharing is caring? Measurement error and the
#' issues arising from combining 3D morphometric datasets. Ecology and
#' Evolution 7:7034-7046.
#'
#' @keywords internal
#' @export
Kmultparallel = function(data, trees, burninpercent = 0, iter = 0,
                         verbose = TRUE) {
    .Defunct(msg =
      "Kmultparallel() is defunct and has been removed from GeometricMorphometricsMix. Please use geomorph::physignal() instead.")
}


#' Print method for parallel_Kmult objects (defunct)
#'
#' @description
#' **Defunct**: Removed from the package along with
#' \code{Kmultparallel()}.
#'
#' @param x Ignored (function is defunct).
#' @param ... Ignored (function is defunct).
#'
#' @return Does not return; always signals a defunct error.
#'
#' @keywords internal
#' @export
print.parallel_Kmult = function(x, ...) {
    .Defunct(msg = paste(
        "print.parallel_Kmult() is defunct and has been removed",
        "from GeometricMorphometricsMix along with Kmultparallel()."
    ))
}


#' Plot method for parallel_Kmult objects (defunct)
#'
#' @description
#' **Defunct**: Removed from the package along with
#' \code{Kmultparallel()}.
#'
#' @param x Ignored (function is defunct).
#' @param alpha Ignored (function is defunct).
#' @param title Ignored (function is defunct).
#' @param x_lab Ignored (function is defunct).
#' @param ... Ignored (function is defunct).
#'
#' @return Does not return; always signals a defunct error.
#'
#' @keywords internal
#' @export
plot.parallel_Kmult = function(x, alpha = 0.25, title = NULL,
                               x_lab = "Kmult", ...) {
    .Defunct(msg = paste(
        "plot.parallel_Kmult() is defunct and has been removed",
        "from GeometricMorphometricsMix along with Kmultparallel()."
    ))
}


#' Summary method for parallel_Kmult objects (defunct)
#'
#' @description
#' **Defunct**: Removed from the package along with
#' \code{Kmultparallel()}.
#'
#' @param object Ignored (function is defunct).
#' @param ... Ignored (function is defunct).
#'
#' @return Does not return; always signals a defunct error.
#'
#' @keywords internal
#' @export
summary.parallel_Kmult = function(object, ...) {
    .Defunct(msg = paste(
        "summary.parallel_Kmult() is defunct and has been removed",
        "from GeometricMorphometricsMix along with Kmultparallel()."
    ))
}
