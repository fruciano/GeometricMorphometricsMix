#' Rescale landmark data based on interlandmark distances
#'
#' Convenience function which rescales a dataset of landmarks
#' based on a vector of distances between two landmarks
#'
#' This function can be useful when one has the distance between two landmarks
#' (e.g., obtained with a caliper), but a scale has not been set
#' when acquiring the data (for instance, a scale bar was missing on photos,
#' so configuration of landmarks are scaled in pixels).
#'
#' @param Data array k x p x n of k landmarks and p dimensions for n observations (specimens)
#' @param lm1,lm2 index of the two landmarks whose inter-landmark distance is known
#' @param lengths vector of n lengths (distances between lm1 and lm2 in the appropriate scale)
#'
#' @return The function outputs an array k x p x n of rescaled landmark coordinates
#'
#'
#' @export
rescale_by_landmark_distance=function(Data, lm1, lm2, lengths) {
  original_distances=unlist(lapply(seq_len(dim(Data)[3]), function(i)
    pdistance(Data[lm1,,i], Data[lm2,,i])
      ))
  # Distances in the data with the original units
  # (e.g., pixels)
  scales=lengths/original_distances
  # Compute scale factors

  Data_rescaled=Data
  for (i in seq_len(dim(Data)[3])) {
    Data_rescaled[,,i]=Data[,,i]*scales[i]
  }
  # rescale data

return(Data_rescaled)
}
