#' Body arching vector from brown trout study
#'
#' A list containing the body arching vector obtained using common principal
#' components following the procedure described in Fruciano et al. 2020, and the
#' GPA consensus used to align individual models when the vector was computed.
#' The GPA consensus is provided to align new data using the same consensus),
#'  but using the consensus for downstream work is not recommended.
#'
#' The object is a list with the following elements:
#' \describe{
#'   \item{GPA_consensus_used}{A matrix of landmark coordinates for the GPA consensus to which the models have been aligned (brown trout dataset in this package)}
#'   \item{arching_vector_CPCA}{A numeric vector of shape change obtained using common principal component analysis as detailed in Fruciano et al. 2020.}
#' }
#'
#' @docType data
#' @usage data(arching_vector)
#' @format A list with two components as described above.
#' @references Fruciano C., Schmidt, I., Ramirez Sanchez, M.M., Morek, W., Avila Valle, Z.A., Talijancic, I., Pecoraro, C., Schermann Legionnet, A. 2020. Tissue preservation can affect geometric morphometric analyses: a case study using fish body shape. Zoological Journal of the Linnean Society 188:148-162.
#' @seealso \code{ProjectOrthogonal} for projection/removal of an effect (see also the functions that implement correction for arching in this package).
#' @keywords datasets
#' @examples
#' data(arching_vector)
#' a = arching_vector
#' names(a)
#' dim(a$GPA_consensus_used)
#' length(a$arching_vector_CPCA)
 #' @name arching_vector
 #' @aliases arching_vector
NULL
