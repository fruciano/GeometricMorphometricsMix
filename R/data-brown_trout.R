#' Brown trout landmark data
#'
#' A list containing landmark coordinates for a subset of brown trout
#' specimens used in Fruciano et al. 2014 (Biological Journal of the
#' Linnean Society) and Fruciano et al. 2020 (Zoological Journal of
#' the Linnean Society). The subset originates from two rivers and was
#' digitised by a single operator.
#'
#' Notice that the dataset has been realigned using Generalized
#' Procrustes Analysis (GPA) and corrected for body arching using the
#' arching vector in \code{data(arching_vector)}.
#' Because this is a small subset of the fish in the original study
#' (Fruciano et al. 2014), and contain only data from a single
#' operator and treatment from Fruciano et al. 2020,
#' the dataset should be strictly considered as a "toy" dataset and
#' not used for formal statistical analyses. For conclusions about
#' biological variation and measurement error due to preservation,
#' please refer to the original studies (Fruciano et al. 2014,
#' Fruciano et al. 2020, respectively).
#'
#' The object is a list with the following components:
#' \describe{
#'   \item{with_arching}{Landmark coordinates for specimens without
#'   correction for body arching.}
#'   \item{without_arching}{Landmark coordinates for the same
#'   specimens after correction for body arching.}
#'   \item{Factors}{A data.frame with metadata for each specimen (e.g.,
#'   centroid size, sex and river).}
#' }
#'
#' @docType data
#' @usage data(brown_trout)
#' @format A list with components described above.
#' @references Fruciano C, Pappalardo AM, Tigano C, Ferrito V. 2014.
#' Phylogeographical relationships of Sicilian brown trout and the
#' effects of genetic introgression on morphospace occupation.
#' Biological Journal of the Linnean Society 112:387-398.
#' @references Fruciano C., Schmidt, I., Ramirez Sanchez, M.M.,
#' Morek, W., Avila Valle, Z.A., Talijancic, I., Pecoraro, C.,
#' Schermann Legionnet, A. 2020. Tissue preservation can affect
#' geometric morphometric analyses: a case study using fish body
#' shape. Zoological Journal of the Linnean Society 188:148-162.
#' @keywords datasets
#' @examples
#' data(brown_trout)
#' d = brown_trout
#' names(d)
#' str(d$Factors)
 #' @name brown_trout
 #' @aliases brown_trout
NULL
