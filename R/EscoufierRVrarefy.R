#' Rarefied version of Escoufier RV coefficient
#'
#' Computes a rarefied estimate of the Escoufier RV coefficient
#' to account for the dependence on sample size
#' of RV, so that RV can be compared among groups with
#' different number of observations (sample sizes)
#'
#' This function computes a rarefied estimate of Escoufier RV coefficient
#' as suggested by Fruciano et al 2013 - Plos One
#' This can be useful to compare RV among groups with the same variables
#' but different sample sizes (as RV depends on sample size, see Fruciano et al 2013,
#' where this procedure is described).
#' The idea is the one rarefies the two groups at the same sample size
#'
#' @section Notice:
#' the function does NOT perform GPA on each rarefied sample
#' this may or may not make a difference in estimates.
#' In many cases, it will probably not make much difference
#' (e.g., Fig. 2 in Fruciano et al 2013)
#'
#' @section Citation:
#' If you use this function please cite both
#' Fruciano et al. 2013 (for using the rarefaction procedure)
#' and Escoufier 1973 (because the procedure is based on Escoufier RV)
#'
#' @param Block1,Block2 Matrices or data frames
#' containing each block of variables
#' (observations in rows, variables in columns).
#' @param reps number of resamplings to obtain the rarefied estimate
#' @param samplesize sample size to which the rarefaction procedure is carried out
#'
#' @seealso \code{\link{EscoufierRV}}
#' @return The function outputs a list with the following elements:
#'  \describe{
#'   \item{Rarefied_RV}{Mean rarefied RV}
#'   \item{Quantiles}{2.5\%, 50\% (median) and 97.5\% percentiles of the rarefied RV}
#'   \item{AllRarefiedSamples}{All RV values obtained using the rarefaction procedure}
#' }
#'
#' @references Escoufier Y. 1973. Le Traitement des Variables Vectorielles. Biometrics 29:751-760.
#' @references Fruciano C, Franchini P, Meyer A. 2013. Resampling-Based Approaches to Study Variation in Morphological Modularity. PLoS ONE 8:e69376.
#' @examples
#' library(MASS)
#' set.seed(123)
#' Pop=mvrnorm(100000,mu=rep(0,100), Sigma=diag(100))
#' # Create a population of 100,000 'individuals'
#' # as multivariate normal random data
#' # We will consider the first 20 columns as the first
#' # block of variables, and the following one as the second block
#'
#' A=Pop[1:50,]
#' B=Pop[501:700,]
#' # Take two groups (A and B)
#' # from the same population (there should be no difference
#' # between them)
#'
#' EscoufierRV(A[,1:20],A[,21:ncol(A)])
#' EscoufierRV(B[,1:20],B[,21:ncol(B)])
#' # Notice how we obtain very different values of Escoufier RV
#' # (this is because they two groups have very different
#' # sample sizes, one 50 observations, the other 200)
#'
#' RarA=RVrarefied(A[,1:20],A[,21:ncol(A)],rep=1000,samplesize=30)
#' RarB=RVrarefied(B[,1:20],B[,21:ncol(A)],rep=1000,samplesize=30)
#' RarA$Rarefied_RV
#' RarB$Rarefied_RV
#' # Rarefying both groups at the same sample size
#' # (in this case 30)
#' # it is clear that the two groups have very similar levels
#' # of association between blocks
#'
#' RarA$Quantiles
#' RarB$Quantiles
#' # And their intevals clearly overlap
#'
#' @export
RVrarefied = function(Block1, Block2, reps = 1000, samplesize) {
         if (nrow(Block1) != nrow(Block2)) {
           stop(paste("Error: the two blocks should have the same number of rows (observations)"))
            }
    BothBlocks = cbind(Block1, Block2)
    endB1 = ncol(Block1)
    startB2 = endB1 + 1
    sizeboth = ncol(BothBlocks)
    RV = vector(length = reps)
    for (i in 1:reps) {
        NewSamp = cbind(Block1, Block2)[sample(1:nrow(Block1), samplesize, replace = TRUE), ]
        RV[i] = EscoufierRV(cbind(NewSamp[, 1:endB1]), cbind(NewSamp[, startB2:sizeboth]))
    }
	
	if(TRUE %in% is.na(RV)) {
	warning(paste("Some of the rarefied RV have given NA \n
		This might be due to random samples all with the same observations\n
		Please, inspect your data and how frequent it is in the RV computed for rarefied samples\n
		to ensure that the results are still meaningful
		"))
	}
	
	RV2=na.omit(RV)
	# remove any NAs
	
    Results = list(Rarefied_RV = mean(RV2), Quantiles = quantile(RV2, c(0.025, 0.5, 0.975)), AllRarefiedSamples = RV)
    return(Results)
}



#' Escoufier RV coefficient
#'
#' Computes the Escoufier RV coefficient
#'
#'
#' This function computes the usual version of the Escoufier RV coefficient (Escoufier, 1973),
#' which quantifies the level of association between two multivariate blocks
#' of variables. The function accepts two blocks of variables, either two data frames
#' or two matrices each of n observations (specimens) as rows.
#' The two blocks must have the same number of rows (specimens), but can have
#' different number of columns (variables, such as landmark coordinates).
#' The Escoufier RV has been shown (Fruciano et al. 2013) to be affected
#' by sample size so comparisons of groups (e.g., species, populations)
#' with different sample size should be avoided, unless steps are taken to account
#' for this problem
#'
#' @param Block1,Block2 Matrices or data frames
#' containing each block of variables
#' (observations in rows, variables in columns).
#'
#' @return The function returns a number, corresponding to the Escoufier RV coefficient
#' @references Escoufier Y. 1973. Le Traitement des Variables Vectorielles. Biometrics 29:751-760.
#' @references Fruciano C, Franchini P, Meyer A. 2013. Resampling-Based Approaches to Study Variation in Morphological Modularity. PLoS ONE 8:e69376.
#'
#' @seealso \code{\link{RVrarefied}}
#'
#' @examples
#' library(MASS)
#' set.seed(123)
#' A=mvrnorm(100,mu=rep(0,100), Sigma=diag(100))
#' # Create a sample of 100 'individuals'
#' # as multivariate normal random data
#' # We will consider the first 20 columns as the first
#' # block of variables, and the following one as the second block
#'
#' EscoufierRV(A[,1:20],A[,21:ncol(A)])
#' # Compute the EscoufierRV using the two blocks of variables
#'
#' @export
EscoufierRV = function(Block1, Block2) {
    if (nrow(Block1) != nrow(Block2)) {
        stop(paste("Error: the two blocks should have the same number of rows (observations)"))
    } else {
        BothBlocks = cbind(Block1, Block2)
        nvarB1 = ncol(Block1)
        COV = cov(BothBlocks)
        RV = sum(diag(COV[1:nvarB1, (nvarB1 + 1):ncol(COV)] %*% COV[(nvarB1 + 1):nrow(COV), 1:nvarB1]))/(sqrt(sum(diag(COV[1:nvarB1,
            1:nvarB1] %*% COV[1:nvarB1, 1:nvarB1])) * sum(diag(COV[(nvarB1 + 1):ncol(COV), (nvarB1 + 1):ncol(COV)] %*% COV[(nvarB1 +
            1):ncol(COV), (nvarB1 + 1):ncol(COV)]))))
        return(RV)
    }
}
