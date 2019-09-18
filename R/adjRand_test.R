#' Test the significance of the adjusted Rand index
#'
#' Permutation test of the adjusted Rand index, which quantifies the level of agreement
#' between two partitions (e.g., two schemes of classification of the same individuals obtained with two methods)
#'
#' The adjusted Rand index (Hubert and Arabie 1985), is an adjusted for chance version of the Rand index (Rand 1971).
#' The adjusted Rand index has an expected value of zero in the case of random partitions,
#' and values approaching one as the two partitions become more similar to each other
#' (with one being perfect match of the classification).
#'  This function implements the permutation test proposed by Qannari et al. (2014)
#' to obtain a p value against the null hypothesis of independence of the two partitions.
#'
#'
#' This function is useful in various contexts, such as in integrative taxonomy
#' when comparing the classification of individual specimens obtained using different data
#' (e.g., sequence data and morphometric data).
#' For an example of the application of this technique with the classification obtained with genetic data and
#' morphometric data for multiple traits, see Fruciano et al. 2016.
#'
#'
#' @section Notice:
#' The function requires internally the package mclust.
#'
#' @param A,B numerical or character vectors reflecting the assignment of individual observations to groups
#' @param perm number of permutations
#' @return The function outputs a vector with the adjusted Rand index and the p value obtained from the permutation test
#'
#' @section Citation:
#' If you use this function in the context of integrative taxonomy or similar
#' (comparison of classification/unsupervised clustering with biological data), please cite all the papers in the references
#' (otherwise, please use the relevant citations for the context).
#'
#' @references Fruciano C, Franchini P, Raffini F, Fan S, Meyer A. 2016. Are sympatrically speciating Midas cichlid fish special? Patterns of morphological and genetic variation in the closely related species Archocentrus centrarchus. Ecology and Evolution 6:4102-4114.
#' @references Hubert L, Arabie P. 1985. Comparing partitions. Journal of Classification 2:193-218.
#' @references Qannari EM, Courcoux P, Faye P. 2014. Significance test of the adjusted Rand index. Application to the free sorting task. Food Quality and Preference 32, Part A:93-97.
#' @references Rand WM. 1971. Objective Criteria for the Evaluation of Clustering Methods. Journal of the American Statistical Association 66:846-850.
#'
#' @examples
#' library(mclust)
#' set.seed(123)
#'
#' irisBIC = mclustBIC(iris[,-5])
#' mclustBIC_classification = summary(irisBIC,iris[,-5])$classification
#' original_classification = iris[,5]
#' # This is one of the examples in the package mclust
#' # Here a classification algorithm is used on the iris dataset
#'
#' adjustedRandIndex(mclustBIC_classification, original_classification)
#' # The mclust package allows computing the adjusted Rand index
#' # which quantifies the agreement between the original (correct) classification
#' # and the one obtained with the algorithm.
#' # However, it is not clear whether the adjusted Rand index is "large enough"
#' # compared to the null hypothesis of independence between the two classification schemes
#'
#' adjRand_test(mclustBIC_classification, original_classification, perm = 999)
#' # For that, we use the function adjRand_test, which performs the permutation test
#' # of Qannari et al. 2014 (in this case p<0.001, as 1000 permutations have been used).
#'
#' adjRand_test(original_classification, original_classification, perm = 999)
#' # As it can be seen, in the ideal case of the exact same grouping,
#' # the adjusted Rand index takes a value of 1 (which is obviously significant)
#'
#'
#' @export
adjRand_test=function(A, B, perm=999) {
  if (length(A)!=length(B)) { stop("A and B should have the same length") }
  # Make sure that the two groups of partitions have the same length

  ARIo=mclust::adjustedRandIndex(A, B)
  # Observed adjusted Rand index

  Aperm=lapply(seq_len(perm), function(X) sample(A, length(A), replace=FALSE))
  Bperm=lapply(seq_len(perm), function(X) sample(B, length(B), replace=FALSE))
  # Generate permuted samples

  ARIperm=unlist(lapply(seq_len(perm),
                        function(i) mclust::adjustedRandIndex(Aperm[[i]], Bperm[[i]])))
  # Compute adjusted Rand index for the permuted samples

  m=mean(ARIperm); v=var(ARIperm)
  # compute mean and variance of the permuted samples

  NARI=(ARIperm-m)/sqrt(v)
  # compute NARI (normalized ARI)

  NARIo=(ARIo-m)/sqrt(v)
  # compute observed NARI

  p_value=(length(which(NARI>NARIo))+1)/(perm+1)
  # Compute p value as proportion of permuted NARI larger than the observed

Results=c(ARIo, p_value)
names(Results)=c("Adjusted_Rand_index", "p_value")
return(Results)
}

