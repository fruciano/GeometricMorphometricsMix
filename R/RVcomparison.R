#' Compare Escoufier RV coefficient between groups
#'
#' Performs permutation tests for the difference in
#' Escoufier RV between groups. For multiple groups, performs
#' pairwise comparisons between all pairs of groups.
#'
#' This function is one of the solutions proposed by Fruciano et al. 2013
#' to deal with the fact that values of Escoufier RV coefficient (Escoufier 1973),
#' which is routinely used to quantify the levels of association between multivariate
#' blocks of variables (landmark coordinates in the case of morphometric data,
#' but it might be any multivariate dataset)
#' are not comparable across groups with different number of observations/individuals
#' (Smilde et al. 2009; Fruciano et al. 2013).
#' The solution is a permutation test. This test was originally implemented by Adriano Franchini
#' in the Java program RVcomparator,
#' of which this function is an updated and improved version
#'
#' @section Notice:
#' The values of RV for each of the groups the function outputs
#' are the observed ones, and can be reported but their value should
#' not be compared. If one wants to obtain comparable values
#' one solution (Fruciano et al 2013) is to obtain rarefied estimates,
#' which can be done with the function RVrarefied in this package
#'
#' @section Citation:
#' If you use this function please cite both
#' Fruciano et al. 2013 (for using the rarefaction procedure)
#' and Escoufier 1973 (because the procedure is based on Escoufier RV)
#'
#' @param Block1,Block2 Matrices or data frames
#' containing each block of variables
#' (observations in rows, variables in columns).
#' @param group factor or vector indicating group membership for observations.
#' @param perm number of permutations for the test
#' @param center whether the groups should be mean-centered prior to the test
#'
#' @seealso \code{\link{EscoufierRV}}
#' @seealso \code{\link{RVrarefied}}
#' @return A data frame with one row per pairwise comparison and the following columns:
#'  \describe{
#'   \item{group1}{Name of the first group in the comparison}
#'   \item{group2}{Name of the second group in the comparison}
#'   \item{Observed_RV_group1}{Observed Escoufier RV for the first group in the comparison}
#'   \item{Observed_RV_group2}{Observed Escoufier RV for the second group in the comparison}
#'   \item{Absolute_difference_in_RV}{Absolute difference in the observed Escoufier RV between the two groups}
#'   \item{p_value}{p value of the permutation test}
#' }
#' For multiple groups, the data frame includes additional columns identifying the groups being compared.
#'
#' @references Escoufier Y. 1973. Le Traitement des Variables Vectorielles. Biometrics 29:751-760.
#' @references Fruciano C, Franchini P, Meyer A. 2013. Resampling-Based Approaches to Study Variation in Morphological Modularity. PLoS ONE 8:e69376.
#' @references Smilde AK, Kiers HA, Bijlsma S, Rubingh CM, van Erk MJ. 2009. Matrix correlations for high-dimensional data: the modified RV-coefficient. Bioinformatics 25:401-405.
#' @examples
#'
#' library(MASS)
#' set.seed(123)
#'
#' # Create sample data with different correlation structures
#' S1 = diag(50)  # Uncorrelated variables for group 1
#' S2 = diag(50)
#' S2[11:50, 11:50] = 0.3  # Some correlation in second block for group 2
#' S2 = (S2 + t(S2)) / 2  # Ensure symmetry
#' diag(S2) = 1
#'
#' # Generate data for two groups
#' A = mvrnorm(30, mu = rep(0, 50), Sigma = S1)
#' B = mvrnorm(40, mu = rep(0, 50), Sigma = S2)
#'
#' # Combine data and create group labels
#' combined_data1 = A[, 1:20]
#' combined_data2 = A[, 21:50]
#' combined_data1 = rbind(combined_data1, B[, 1:20])
#' combined_data2 = rbind(combined_data2, B[, 21:50])
#' groups = c(rep("GroupA", 30), rep("GroupB", 40))
#'
#' # Perform RV comparison
#' result = RVcomparison(combined_data1, combined_data2, group = groups, perm = 99)
#' print(result)
#'
#' # Example with three groups for pairwise comparisons
#' C = mvrnorm(25, mu = rep(0, 50), Sigma = diag(50))
#' combined_data1_3grp = rbind(combined_data1, C[, 1:20])
#' combined_data2_3grp = rbind(combined_data2, C[, 21:50])
#' groups_3 = c(groups, rep("GroupC", 25))
#'
#' result_3grp = RVcomparison(combined_data1_3grp, combined_data2_3grp, 
#'                           group = groups_3, perm = 99)
#' print(result_3grp)
#'
#' @import stats
#' @export
RVcomparison = function(Block1, Block2, group, perm = 999, center = TRUE) {

  # Validate inputs
  validated = validate_RVcomparison_inputs(Block1, Block2, group, perm)
  Block1 = validated$Block1
  Block2 = validated$Block2
  group = validated$group
  n_obs = validated$n_obs

  
  unique_groups = levels(group)
  n_groups = length(unique_groups)
  
  # For single group, return message indicating no comparison possible
  if (n_groups == 1) {
    stop("At least two groups are required for RV comparison")
  }
  
  # Generate all pairwise combinations
  group_pairs = combn(unique_groups, 2, simplify = FALSE)
  n_comparisons = length(group_pairs)
  
  # Initialize result data frame
  results = data.frame(
    group1 = character(n_comparisons),
    group2 = character(n_comparisons),
    Observed_RV_group1 = numeric(n_comparisons),
    Observed_RV_group2 = numeric(n_comparisons),
    Absolute_difference_in_RV = numeric(n_comparisons),
    p_value = numeric(n_comparisons),
    stringsAsFactors = FALSE
  )
  
  # Perform pairwise comparisons
  for (i in seq_len(n_comparisons)) {
    grp1_name = group_pairs[[i]][1]
    grp2_name = group_pairs[[i]][2]
    
    # Extract data for the two groups
    grp1_idx = which(group == grp1_name)
    grp2_idx = which(group == grp2_name)
    
    Group1Block1 = Block1[grp1_idx, , drop = FALSE]
    Group1Block2 = Block2[grp1_idx, , drop = FALSE]
    Group2Block1 = Block1[grp2_idx, , drop = FALSE]
    Group2Block2 = Block2[grp2_idx, , drop = FALSE]
    
    ngrp1 = nrow(Group1Block1)
    ngrp2 = nrow(Group2Block1)
    nvarblk1 = ncol(Group1Block1)
    nvarblk2 = ncol(Group1Block2)
    obstot = ngrp1 + ngrp2
    vartot = nvarblk1 + nvarblk2
    
    # Center data if requested
    if (center == TRUE) {
      Group1 = scale(cbind(Group1Block1, Group1Block2), center = TRUE, scale = FALSE)
      Group2 = scale(cbind(Group2Block1, Group2Block2), center = TRUE, scale = FALSE)
      Both = rbind(Group1, Group2)
      Group1Block1 = Group1[, 1:nvarblk1, drop = FALSE]
      Group1Block2 = Group1[, (nvarblk1 + 1):vartot, drop = FALSE]
      Group2Block1 = Group2[, 1:nvarblk1, drop = FALSE]
      Group2Block2 = Group2[, (nvarblk1 + 1):vartot, drop = FALSE]
    } else {
      Both = cbind(rbind(Group1Block1, Group2Block1), rbind(Group1Block2, Group2Block2))
    }
    
    # Calculate observed RV values
    ObsRVgrp1 = EscoufierRV(Group1Block1, Group1Block2)
    ObsRVgrp2 = EscoufierRV(Group2Block1, Group2Block2)
    obsAbsDiffRV = abs(ObsRVgrp1 - ObsRVgrp2)
    
    # Perform permutation test
    PermutedSets = lapply(seq_len(perm), function(x)
      Both[sample(seq_len(obstot), size = obstot, replace = FALSE), ]
    )
    
    RVdiffperm = unlist(lapply(PermutedSets, function(x)
      abs(
        EscoufierRV(x[1:ngrp1, 1:nvarblk1, drop = FALSE], x[1:ngrp1, (nvarblk1 + 1):vartot, drop = FALSE]) -
          EscoufierRV(x[(ngrp1 + 1):obstot, 1:nvarblk1, drop = FALSE], x[(ngrp1 + 1):obstot, (nvarblk1 + 1):vartot, drop = FALSE])
      )
    ))
    
    pvalue = (length(which(RVdiffperm >= obsAbsDiffRV)) + 1) / (perm + 1)
    
    # Store results
    results$group1[i] = grp1_name
    results$group2[i] = grp2_name
    results$Observed_RV_group1[i] = ObsRVgrp1
    results$Observed_RV_group2[i] = ObsRVgrp2
    results$Absolute_difference_in_RV[i] = obsAbsDiffRV
    results$p_value[i] = pvalue
  }
  
  return(results)
}

# Internal validation function for RVcomparison inputs
validate_RVcomparison_inputs = function(Block1, Block2, group, perm) {
  # Basic dimensional validation
  if (nrow(Block1) != nrow(Block2)) {
    stop("Block1 and Block2 must have the same number of observations (rows)")
  }
  
  # Validate perm parameter
  if (!is.numeric(perm) || length(perm) != 1 || perm != round(perm) || perm <= 0) {
    stop("perm must be a positive integer")
  }
  
  # Convert to data frames if they aren't already
  Block1 = as.data.frame(Block1)
  Block2 = as.data.frame(Block2)
  
  # Check for NA values
  if (any(is.na(Block1))) {
    stop("Block1 contains NA values")
  }
  if (any(is.na(Block2))) {
    stop("Block2 contains NA values")
  }
  
  # Check that blocks contain only numeric data
  if (!all(vapply(Block1, is.numeric, FUN.VALUE = logical(1)))) {
    stop("Block1 must contain only numeric variables")
  }
  if (!all(vapply(Block2, is.numeric, FUN.VALUE = logical(1)))) {
    stop("Block2 must contain only numeric variables")
  }
  
  n_obs = nrow(Block1)
  
  # Validate group input
  if (length(group) != n_obs) {
    stop("Length of group vector must equal number of observations")
  }
  if (any(is.na(group))) {
    stop("group vector contains NA values")
  }
  
  # Return validated data
  return(list(
    Block1 = Block1,
    Block2 = Block2,
    group = as.factor(group),
    n_obs = n_obs
  ))
}

