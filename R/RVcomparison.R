#' Compare Escoufier RV coefficient between two groups
#'
#' Performs a permutation test for the difference in
#' Escoufier RV between two groups
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
#' @param Group1Block1,Group1Block2,Group2Block1,Group2Block2 Matrices or data frames
#' containing the data for each group of observations (specimens) and each block of variables
#' (observations in rows, variables in columns).
#' @param perm number of permutations for the test
#' @param center whether the two groups should be mean-centered prior to the test
#'
#' @seealso \code{\link{EscoufierRV}}
#' @seealso \code{\link{RVrarefied}}
#' @return The function outputs a list with the following elements:
#'  \describe{
#'   \item{ObservedRVgrp1}{Observed Escoufier RV for group 1}
#'   \item{ObservedRVgrp2}{Observed Escoufier RV for group 2}
#'   \item{ObservedAbsoluteDifference}{Absolute difference in the observed Escoufier RV between the two groups}
#'   \item{pvalue}{p value of the permutation test}
#' }
#'
#' @references Escoufier Y. 1973. Le Traitement des Variables Vectorielles. Biometrics 29:751-760.
#' @references Fruciano C, Franchini P, Meyer A. 2013. Resampling-Based Approaches to Study Variation in Morphological Modularity. PLoS ONE 8:e69376.
#' @references Smilde AK, Kiers HA, Bijlsma S, Rubingh CM, van Erk MJ. 2009. Matrix correlations for high-dimensional data: the modified RV-coefficient. Bioinformatics 25:401-405.
#' @examples
#' library(MASS)
#' if(!require(lmf)){ install.packages("lmf") }
#' library(lmf)
#'
#' set.seed(123)
#'
#' S1=diag(100)
#' # Use uncorrelated variables for the first group
#' 
#' S2=diag(100)
#' S2[21:100,21:100]=rnorm((100-20)^2, 0.1, 0.03^2)
#' S2[1:20,1:20]=rnorm(20^2,0.01,0.03^2)
#' lower=S2[lower.tri(S2)]
#' S2=t(S2)
#' S2[lower.tri(S2)]=lower
#' diag(S2)=1
#' S2=nearPD(S2)
#' # Create an ad hoc covariance matrix for the second group
#' 
#' 
#' A=mvrnorm(50, mu=rep(0,100), Sigma=S1)
#' B=mvrnorm(100, mu=rep(10,100), Sigma=S2)
#' # Create two multivariate normal random datasets
#' # using the covariance matrices above
#' # We are going to use the first 20 variables as block 1
#' # and the subsequent variables as block 2
#' 
#' EscoufierRV(A[,1:20], A[,21:ncol(A)])
#' EscoufierRV(B[,1:20], B[,21:ncol(B)])
#' # Notice that the two observed Escoufier RV are different
#' # However, this difference might be due to the
#' # different number of observations for the two groups
#' 
#' RVcomparison(A[,1:20], A[,21:ncol(A)], B[,1:20], B[,21:ncol(B)])
#' # One of the solutions is to perform a permutation test,
#' # which shows us that these differences are actually significant
#' # (another, complementary, solution is to compute rarefied estimates
#' # using the function EscoufierRVrarefy)
#' RVrarefied(A[,1:20], A[,21:ncol(A)], samplesize=40)$Rarefied_RV
#' RVrarefied(B[,1:20], B[,21:ncol(B)], samplesize=40)$Rarefied_RV
#'
#' @export
RVcomparison=function(Group1Block1,Group1Block2,Group2Block1,Group2Block2, perm=999, center=TRUE) {
  
  ngrp1=nrow(Group1Block1)
  ngrp2=nrow(Group2Block1)
  nvarblk1=ncol(Group1Block1)
  nvarblk2=ncol(Group1Block2)
  obstot=ngrp1+ngrp2
  vartot=nvarblk1+nvarblk2
  # Compute a few useful quantities
  
  if (nrow(Group1Block1)!=nrow(Group1Block2)) {
    stop("The first group has different number of observations (specimens) across the two blocks")
    } 
  if (nrow(Group2Block1)!=nrow(Group2Block2)) {
    stop("The second group has different number of observations (specimens) across the two blocks")
    }
  if (ncol(Group2Block1)!=nvarblk1) {
    stop("The two groups must have the same number of variables in each block")
  }
  if (ncol(Group2Block2)!=nvarblk2) {
    stop("The two groups must have the same number of variables in each block")
  }
 #  Checks that the data is consistent across blocks/groups
  
  if (center==TRUE) {
    Group1=scale(cbind(Group1Block1, Group1Block2), center=TRUE, scale=FALSE)
    Group2=scale(cbind(Group2Block1, Group2Block2), center=TRUE, scale=FALSE)
    Both=rbind(Group1, Group2)
    Group1Block1=Group1[,1:nvarblk1]
    Group1Block2=Group1[,(nvarblk1+1):vartot]
    Group2Block1=Group2[,1:nvarblk1]
    Group2Block2=Group2[,(nvarblk1+1):vartot]
    # Scale and prepare data    
    } else {
      Both=cbind(rbind(Group1Block1,Group2Block1),rbind(Group1Block2,Group2Block2))
    # Combine data and compute basic info
    }

  ObsRVgrp1=EscoufierRV(Group1Block1, Group1Block2)
  ObsRVgrp2=EscoufierRV(Group2Block1, Group2Block2)
  obsAbsDiffRV=abs(ObsRVgrp1-ObsRVgrp2)    
  # Observed RV and their difference

  PermutedSets=lapply(seq_len(perm), function(x)
      Both[sample(seq_len(obstot), size=obstot, replace=FALSE),]
          )
  # Create permuted datasets
  RVdiffperm=unlist(lapply(PermutedSets, function(x)
              abs(
                EscoufierRV(x[1:ngrp1,1:nvarblk1],x[1:ngrp1,(nvarblk1+1):vartot])-
                  EscoufierRV(x[(ngrp1+1):obstot,1:nvarblk1],x[(ngrp1+1):obstot,(nvarblk1+1):vartot])  
                )              
              ))

  pvalue=length(which(RVdiffperm>obsAbsDiffRV))/(perm+1)

return(list(ObservedRVgrp1=ObsRVgrp1,
            ObservedRVgrp2=ObsRVgrp2,
            ObservedAbsoluteDifference=obsAbsDiffRV,
            pvalue=pvalue
            ))
}

