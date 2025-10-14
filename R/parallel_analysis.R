#' Perform parallel analysis
#'
#' Parallel analysis based on permutations
#'
#' The function allows performing parallel analysis, which is a way to test for
#' the number of significant eigenvalues/axes in a PCA.
#' In this implementation, a null distribution of eigenvalues is obtained
#' by randomly permuting observations independently for each of the starting variables.
#' To compute p values, the observed eigenvalues are
#' compared to the corresponding eigenvalues from this null distribution.
#'
#' Parallel analysis may be used for dimensionality reduction, retaining
#' only the first block of consecutive significant axes.
#' That is, if for example the first 3 axes were significant, then the fourth not significant,
#' one would keep only the first 3 axes (regardless of significance of the axes from the fifth on).
#' Similarly, if the first axis is not significant, this may suggest lack of a clear structure in the data.
#'
#'
#'
#' The function internally employs three possible strategies to obtain eigenvalues (argument of fun):
#'
#' \itemize{
#'   \item "prcomp" - the function prcomp (default)
#'   \item "fastSVD" - an approach based on the function fast.svd (requires the package corpcor)
#'   \item "shrink" - a decomposition of the covariance matrix estimated using linear shrinkage (much slower, requires the package nlshrink; Ledoit & Wolf 2004)
#' }
#'
#' This choice should not make much difference in terms of the final result.
#' However, for consistency, it is a good idea to use for parallel analysis the same function
#' used for the actual PCA (this is why these three options are provided).
#'
#'
#'
#'
#' @param X Matrix or data frame containing the original data
#' (observations in rows, variables in columns).
#' @param perm number of permutations
#' @param fun function to use internally to obtain eigenvalues (see Details)
#'
#'
#' @return Vector of class "parallel_analysis" containing the p values
#' for each of the axes of a PCA on the data provided
#'
#'
#' The object of class parallel_analysis returned by the function has a summary() method associated to it.
#' This means that using summary() on an object created by this function,
#' a suggestion on the number of significant axes (if any) is provided (see examples).
#'
#' @section Citation:
#' The most appropriate citation for this approach to parallel analysis (using permutations) is Buja & Eyuboglu (1992).
#'
#'
#' @references Ledoit O, Wolf M. 2004. A well-conditioned estimator for large-dimensional covariance matrices. Journal of Multivariate Analysis 88:365-411.
#' @references Buja A, Eyuboglu N. 1992. Remarks on parallel analysis. Multivariate Behavioral Research 27(4):509-540.
#'
#'
#' @examples
#' set.seed(666)
#' X=MASS::mvrnorm(100, mu=rep(0, 50), Sigma=diag(50))
#' # Simulate a multivariate random normal dataset
#' # with 100 observations and 50 indipendent variables
#'
#' PA=parallel_analysis(X, perm = 999, fun = "fastSVD")
#' # Perform parallel analysis
#'
#' summary(PA)
#' # Look at a summary of the results from parallel analysis
#' # Notice that no axis is significant
#' # This is correct, as we had simulated data with no structure
#'
#' print(PA)
#' # Look at the p values for each axis
#'
#'
#' @export

parallel_analysis=function(X, perm=999, fun=c("prcomp", "fastSVD", "shrink")) {

  fun=ifelse(fun %in% c(c("prcomp", "fastSVD", "shrink")), fun[1], "prcomp")[1]
  if (length(intersect(class(X), c("data.frame", "matrix")))==0) {
    stop("X should be either a matrix or a data.frame")
    }
  if (any(dim(X)<=1)) {
    stop("X should have multiple rows and columns")
  }
  # Initial checks

  Results=pa_PCA(X, perm=perm, method=fun)

return(Results)
}

#' @export
summary.parallel_analysis = function (object, ...) {
  cat("Parallel analysis")
  cat("\n")
  cat("======================")
  cat("\n")
  cat(paste0("Number of tested axes: ", length(object)))
  cat("\n")
   sugg_axes=NA
   if (all(object>0.05)) {sugg_axes=0}
   else if (all(object<=0.05)) {sugg_axes=length(object)}
   else {sugg_axes=(which(object>0.05)[1])-1}
   cat(paste0("Number of axes significant at the 5% level: ", sugg_axes))
   cat("\n")
   cat("======================")
}


# Create dataset for Parallel Analysis (permutation of variable rows
# separate for each variable)
# (not exported)
pa_perm_scheme=function(X, perm=999) {
  nobs=nrow(X)
  PA_PermX=lapply(seq(perm), function(rp) {
    Xp=X
    for (cl in seq_len(ncol(X))) {
      Xp[,cl]=X[sample(seq(nobs), size = nobs), cl]
    }
    return(Xp)
  })
}


# Function to perform parallel analysis in PCA
# (not exported)
pa_PCA=function(X, perm=99, method="prcomp") {
  PA_perm_X=pa_perm_scheme(X, perm = perm)
  # Get permuted samples for parallel analysis

  if (method=="prcomp") {
    obs_PCA=prcomp(X)
    obs_eival=obs_PCA$sdev^2
    PA_perm_eigval=do.call("cbind",
                           lapply(PA_perm_X, function(Y)
                             prcomp(Y)$sdev^2)
    )
  } else if (method=="shrink") {
    obs_PCA=PCA_shrink(X)
    obs_eival=obs_PCA$Eigenvalues
    PA_perm_eigval=do.call("cbind",
                           lapply(PA_perm_X, function(Y)
                             PCA_shrink(Y)$Eigenvalues)
    )
  } else if (method=="fastSVD") {
    obs_PCA=PCA_fastsvd(X)
    obs_eival=obs_PCA$Eigenvalues
    PA_perm_eigval=do.call("cbind",
                           lapply(PA_perm_X, function(Y)
                             PCA_fastsvd(Y)$Eigenvalues)
    )
  }
  # Get observed PCA and eigenvalues for the permuted samples
  # using one of three different methods

  P_values=unlist(lapply(seq_along(obs_eival), function(eiv) {
    (length(which(PA_perm_eigval[eiv,]>=obs_eival[eiv]))+1)/(perm+1)
  }
  ))
  # Compute p values
class(P_values)=c("parallel_analysis", "numeric")
return(P_values)
}



# PCA with fast.svd function
# (not exported)
PCA_fastsvd=function(X) {
  X=scale(X, center = TRUE, scale=FALSE)
  SVD=corpcor::fast.svd(X)
  Eigenvalues=(SVD$d^2)/(nrow(X)-1)
  Vectors=SVD$v
  Scores=X%*%SVD$v

  Results=list(Scores=Scores,
               Vectors=Vectors,
               Eigenvalues=Eigenvalues)
}

# PCA using linear shrinkage estimator of the covariance matrix
# (not exported)
PCA_shrink=function(X) {
  X=scale(X, center=TRUE, scale=FALSE)
  covXshrink=nlshrink::linshrink_cov(X)
  eig_covX=eigen(covXshrink)
  Results=list(Scores=X%*%eig_covX$vectors,
               Eigenvalues=eig_covX$values,
               Eigenvectors=eig_covX$vectors)
}


