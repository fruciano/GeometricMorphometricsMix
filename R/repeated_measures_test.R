#' Perform test on two repeated measures
#'
#' Test based on Hotelling's T squared for the null hypothesis of no effect
#' between two repeated measures (e.g., treatment/control)
#'
#' The function assumes that each individual observation (e.g., specimen) has been measured two times
#' (e.g., at two time points, or between two treatments).
#'
#' If rnames is TRUE (default), the rownames of the matrix or the names along
#' the 3rd dimension (for arrays) will be used to match the order of observations (e.g., specimens)
#' between the two datasets. Otherwise, the function will assume that the observations in T1 and T2
#' are in the same order.
#'
#' This function is useful in various contexts, such as:
#'  \itemize{
#'   \item testing the effect of preservation (Fruciano et al. - under review)
#'   \item testing for variation through time
#' }
#'
#' For instance, in the context of the effect of preservation on geometric morphometrics,
#' it has been argued (Fruciano, 2016) that various studies have improperly used on repeated measures data
#' methods developed for independent observations, and this can lead to incorrect inference.
#'
#' @section Notice:
#' The function requires internally non-singular matrices
#' (for instance, number of cases should be larger than the number of variables).
#' One solution can be to perform a principal component analysis and use the scores
#' for all the axes with non-zero and non-near-zero eigenvalues.
#' To overcome some situations where a singular matrix can occurr, the function can
#' use internally a shrinkage estimator of the covariance matrix (Ledoit & Wolf 2004).
#' This is called setting shrink = TRUE.
#' However, in this case, the package nlshrink should have been installed.
#'
#' @param T1,T2 matrices (n x p of n observation for p variables)
#' or arrays (t x p x n of n observations, t landmarks in p dimensions),
#' @param rnames if TRUE (default) the rownames of the matrix or the names along
#' the 3rd dimension (for arrays) will be used to match the order
#' @param shrink if TRUE, a shrinkage estimator of covariance is used internally
#' @return The function outputs a matrix n x p of the original data projected
#' to the subspace orthogonal to the vector
#'
#' @references Fruciano C. 2016. Measurement error in geometric morphometrics. Development Genes and Evolution 226:139-158.
#' @references Fruciano C., Schmidt, I., Ramirez Sanchez, M.M., Morek, W., Avila Valle, Z.A., Talijancic, I., Pecoraro, C., Schermann Legionnet, A. under review.
#' Geometric morphometric analyses can be affected by tissue preservation: a case study using fish body shape.
#' @references Ledoit O, Wolf M. 2004. A well-conditioned estimator for large-dimensional covariance matrices. Journal of Multivariate Analysis 88:365-411.
#'
#' @examples
#' library(MASS)
#' set.seed(123)
#' A=mvrnorm(50,mu=rep(1,10),Sigma=diag(10))
#' # Generate a random dataset from a multivariate normal distribution
#' # (50 observations, 10 variables)
#' NoEff=do.call("cbind", lapply(1:ncol(A), function(x)
#'                        rnorm(nrow(A), mean=0, sd=0.1)) )
#' # Generate some random noise at each variable
#'
#' B=A+NoEff
#' # Create a new matrix by summing the original matrix to the random noise
#' # (simulating the case of no effect)
#'
#' repeated_measures_test(T1=A, T2=B, rnames=FALSE)
#' # The comparison is not significant
#' # (this is the expected result, as we have simulated
#' # the case of no effect between the two repeated measures)
#' @export
repeated_measures_test=function(T1, T2, rnames=TRUE, shrink=FALSE) {
	if (class(T1)=="array") {
		T1=Morpho::vecx(T1)
			}
	if (class(T2)=="array") {
		T2=Morpho::vecx(T2)
			}
	# If an array, transform in a 2D matrix
	if (rnames==TRUE) {
	T2=T2[rownames(T1),]
	warning("The names of the observations in the two datasets will be used for matching them")
	} else {
	warning("Names are not used, the observations will be assumed to be in the same order")
	}
	# If rnames is TRUE, reorder the second array
	# based on the rownames of the first array
	if (nrow(T1)!=nrow(T2)) {
	stop(paste("The two sets have different number of rows (observations)"))
	}
	if (ncol(T1)!=ncol(T2)) {
	stop(paste("The two sets have different number of columns (variables)"))
	}
	if (nrow(T1)<=ncol(T1)) {
			warning('Number of cases less or equal to the number of variables')
			}
	ObsEuclideandD=dist(rbind(colMeans(T1),colMeans(T2)))
	ObsPairedT2=HotellingT2p(T1,T2, shrink=shrink)
	# Compute  Euclidean distances between the two treatments and
	# paired Hotteling T squared from one observation
	# in one treatment and the same observation in the second treatment
	Results=c(ObsEuclideandD,
			ObsPairedT2$HottelingT2,
			ObsPairedT2$Fstat,
			ObsPairedT2$p_value
			)
	names(Results)=c("EuclideanD",
					"HotellingT2",
					"Fstat",
					"p_value"
				)
return(Results)
}


# Function to resample within subjects
# (not exported)
resamplewithin2=function(A1, A2) {
id1=sample(1:2, nrow(A1), replace=TRUE)
id2=sapply(id1, function(x) if (x==1) 2 else 1)
Alist=list(A1,A2)
A1res=do.call("rbind", lapply(seq_len(length(id1)), function(x) Alist[[id1[x]]][x,]))
A2res=do.call("rbind", lapply(seq_len(length(id2)), function(x) Alist[[id2[x]]][x,]))
Ares=list(A1res,A2res)
return(Ares)
}


# Function to compute Hotteling T squared
# with repeated measures data
# (not exported)
HotellingT2p=function(A1, A2, shrink=FALSE) {
	D=A2-A1
	Di=colMeans(D)
	if (shrink==TRUE) {
	  S=nlshrink::linshrink_cov(D)
	} else {
	  S=cov(D)
	}
	n=nrow(A1)
	p=ncol(A1)
	df2=n-p
	pT2=n*rbind(Di)%*%solve(S)%*%cbind(Di)
		if (n<=p) {
			Fstat=NA
			pval=NA
			# warning('More variables than cases,
			# no F statistic nor p-value will be computed')
			} else {
	Fstat = pT2 / (p * (n-1) / df2)
	pval = 1 - pf(Fstat, df1=p, df2=df2)
		}
Results=list(HottelingT2=pT2,
				Fstat=Fstat,
				p_value=pval)

return(Results)
}

