#' Check the relative positions for a set of landmarks, compared to a reference specimen
#'
#' The function compares the relative position of a set of landmarks to the one observed in a reference specimen.
#' It also outputs which coordinates do not much the pattern observed in the reference specimen.
#'
#' Compare relative positions of landmarks to the one observed in a reference specimen.
#' This function is useful to identify specimens with switched landmark positions and similar problems
#' when a dataset has been collected using consistent criteria
#' (e.g., a set of fish body pictures, with the mouth-tail axis approximately horizontal for all specimens).
#'
#' For instance, if we want to check that landmarks 1, 2, and 3 are all aligned along the y coordinate in 2D,
#' we will use a specimen (usually the first), checking that it has been landmarked correctly.
#' Then, we will run the function on landmarks 1, 2, and 3 (by setting LM_to_check=c(1, 2, 3)).
#' The function will output a list of which specimens seem to be problematic and which landmarks and coordinates for these specimens seem problematic.
#' The putatively problematic coordinates will be indicated with "0".
#' Clearly, if we are only interested in 1, 2, and 3 being along y, it will not matter if in some case
#' we will find "0" under Coord1, but only if we find a "0" under Coord2.
#'
#' @section Notice:
#' Generally, it is better to use this function with a small number
#' of carefully selected landmarks, instead of running it on all landmarks at the same time.
#'
#' The parameter only_by_landmark_order allows to set whether all
#' combinations of landmarks should be tested (default), or not.
#'
#' If only_by_landmark_order is set to TRUE, each landmark in LM_to_check
#' will be only tested with the next one. For instance if
#' LM_to_check=c(1,2,3) the function will only compare the first landmark with the second,
#' and the second with the third.
#'
#' If the function cannot find potentially problematic specimens, a message
#' to this effect will be presented.
#'
#' @param Dataset k x p x n array
#' of k landmarks and p dimensions for n observations (specimens).
#' @param LM_to_check Vector with the indices of which landmarks of Dataset should be tested
#' @param reference_specimen The index of the specimen to use as reference (by default, the first).
#' @param only_by_landmark_order If TRUE, each landmark in LM_to_check will be tested relative to the next in LM_to_check;
#' otherwise all the pairwise relative positions between landmarks will be run
#' @param use_specimen_names whether the names of specimens in Dataset should be used for the output
#'
#'
#' @return The function outputs a list with the following elements:
#'  \describe{
#'   \item{potentially_problematic_idx}{Indices of potentially problematic specimens}
#'   \item{potentially_problematic_LMs}{List of potentially problematic landmarks and coordinates
#'   for the specimens in potentially_problematic_idx}
#' }
#' @importFrom utils combn
#' @export
LM_relativepos_check=function(Dataset, LM_to_check,
                              reference_specimen=1,
                              only_by_landmark_order=FALSE,
                              use_specimen_names=TRUE) {
  if (is.array(Dataset)==FALSE) stop("Dataset should be a k x p x n array of k landmarks and p dimensions for n observations (specimens)")
  if (!dim(Dataset)[2] %in% c(2, 3)) stop("The function accepts only 2D and 3D data")


  if (only_by_landmark_order==TRUE) {
    LM_combinations=cbind(LM_to_check[seq_len(length(LM_to_check)-1)],
                          LM_to_check[2:length(LM_to_check)])
  } else {
    LM_combinations=t(combn(LM_to_check, 2))
  }
  colnames(LM_combinations)=c("LM1", "LM2")
  # Get combinations of landmarks


  relative_positions=array(NA, dim = c(dim(LM_combinations), dim(Dataset)[3]))
  for (i in seq(dim(Dataset)[3])) {
    relative_positions[,,i]=combinations_higher(comb = LM_combinations,
                                                specimen = Dataset[,,i])
  }
  # For each specimen, make a comparison for each landmark combinations/coordinate
  # (for the first being larger than the second)

  evaluation=lapply(seq(dim(relative_positions)[3]), function(i)
    relative_positions[,,i]==relative_positions[,,reference_specimen])

  if (all(unlist(lapply(evaluation, all)))) {
    message("Cannot find potentially problematic specimens")
  } else {
    potentially_problematic_idx=which(unlist(lapply(evaluation, all))==FALSE)
    potentially_problematic_LMs=lapply(evaluation[potentially_problematic_idx],
                function(X) {
                  if (length(LM_to_check)==2) {
                  T1=c(LM_to_check, X)
                  names(T1)=c("LM1", "LM2",
                              paste0("Coord", seq_len(dim(Dataset)[2])))
                  } else {
      T1=cbind(rbind(LM_combinations[apply(X, 1, function(rw)
                                            any(rw==FALSE)),]),
            rbind(X[apply(X, 1, function(rw) any(rw==FALSE)),]))
      colnames(T1)[3:ncol(T1)]=paste0("Coord", seq_len(ncol(T1)-2))
                  }
      return(T1)
                }
      )

    if (use_specimen_names==TRUE) {
      names(potentially_problematic_LMs)=dimnames(Dataset)[[3]][potentially_problematic_idx]
    } else {
      names(potentially_problematic_LMs)=potentially_problematic_idx
    }
    return(list(potentially_problematic_idx=potentially_problematic_idx,
                potentially_problematic_LMs=potentially_problematic_LMs))
  }
}

combinations_higher=function(comb, specimen) {
  out=matrix(NA, nrow(comb), dim(specimen)[2])
  for (i in seq_len(nrow(comb))) {
    out[i,]=specimen[comb[i,1],]> specimen[comb[i,2],]
  }
  return(out)
}
