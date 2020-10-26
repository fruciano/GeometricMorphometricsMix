#' Compute scaled variance of eigenvalues
#'
#' Compute estimates of the scaled variance of eigenvalues using only the positive eigenvalues
#'
#' The function allows computing the scaled variance of eigenvalues (Pavlicev et al. 2009)
#' under a variety of settings.
#' The scaled variance of eigenvalues is a commonly used index of morphological integration.
#' Its value is comprised between 0 and 1, with larger values suggesting stronger integration.
#'
#' Only positive eigenvalues are used in the computations used in this function.
#'
#'
#' The function employes two possible strategies to obtain eigenvalues:
#'
#' \itemize{
#'   \item a singular value decomposition of the data matrix (default)
#'   \item an eigenvalue decomposition of the covariance matrix estimated using linear shrinkage (option shrinkage=TRUE; Ledoit & Wolf 2004)
#' }
#'
#' Further, the function allows obtaining bootstrapped estimates by
#' setting boot to the number of bootstrap replicates (resampling with replacement)
#'
#' It is also possible to obtain rarefied estimates by setting rarefy to the desired sample size.
#' This is useful when comparing the scaled variance of eigenvalues across multiple groups with different sample sizes.
#' In this case, the suggestion is to use either the smallest sample size or less
#'
#' Using a bootstrap estimate with the singular value decomposition approach
#' represents a good compromise between computation time and accuracy.
#' This should be complemented by rarefaction to the smallest sample size (or lower)
#' in case one wants to compare the value obtained across different groups.
#'
#' @section Notice:
#' When boot>0 the rarefied estimates are based on sampling with replacement (bootstrap).
#' However, if boot=0, then a \strong{single} rarefied estimate is obtained by sampling without replacement.
#' In this case, the user should repeat the same operation multiple times (e.g., 100)
#' and take the average of the scaled variance of eigenvalues obtained.
#'
#' Also notice that using the shrinkage-based estimation of the covariance matrix requires longer computational time and memory.
#' This option requires the package \emph{nlshrink}
#'
#' @section Computational details:
#' For the computation, this function uses only positive eigenvalues (which are also used to identify data dimensionality).
#' The eigenvalues are first scaled by dividing them by their sum (Young 2006),
#' then their variance is computed as population variance (rather than sample variance; see Haber 2011).
#' Finally, this value is standardized to a scale between 0 and 1 by dividing it by its maximum theoretical value
#' of (p-1)/p^2 (where p is the number of dimensions) - this is the same scaling used in the software MorphoJ (Klingenberg 2011).
#'
#'
#' @param data_matrix Matrix or data frame containing the original data
#' (observations in rows, variables in columns).
#' @param boot number of bootstrap resamples (no bootstrap if 0)
#' @param rarefy either a logical to determine whether rarefaction will be performed
#' or a number indicating the sample size at which the samples will be rarefied
#' @param shrinkage logical on whether the analysis should be based on a covariance matrix obtained through linear shrinkage
#'
#'
#' @return If boot=0 the function outputs a vector containing
#' the scaled variance of eigenvalues
#' and the number of dimensions used in the computations.
#' If, instead, boot>0 (recommended) the function outputs a list containing
#' \itemize{
#'   \item the mean scaled variance of eigenvalues across all bootstrap samples
#'   \item the median number of dimensions used across all bootstrap samples
#'   \item the 95% confidence intervals for scaled variance of eigenvalues and dimensions
#'   \item the scaled variance of eigenvalues and dimensions used for each of the bootstrap replicates
#' }
#'
#' @references Ledoit O, Wolf M. 2004. A well-conditioned estimator for large-dimensional covariance matrices. Journal of Multivariate Analysis 88:365-411.
#' @references Young NM. 2006. Function, ontogeny and canalization of shape variance in the primate scapula. Journal of Anatomy 209:623-636.
#' @references Pavlicev M, Cheverud JM, Wagner GP. Measuring Morphological Integration Using Eigenvalue Variance. Evolutionary Biology 36(1):157–170.
#' @references Haber A. 2011. A Comparative Analysis of Integration Indices. Evolutionary Biology 38:476-488.
#' @references Klingenberg CP. 2011. MorphoJ: an integrated software package for geometric morphometrics. Molecolar Ecology Resources 11:353-357.
#'
#' @importFrom corpcor fast.svd
#' @import stats
#'
#' @export

scaled_variance_of_eigenvalues=function(data_matrix, boot=999,
                                        rarefy=FALSE, shrinkage=FALSE) {

	if (class(data_matrix) %in% c("data.frame", "matrix")==FALSE) {stop("data_matrix should be either a matrix or a data frame")}

  # The first option is the normal scaled variance of eigenvalues,
  # but based only on the positive eigenvalues obtained from corpcor::fast.svd
  if (all(c(rarefy==FALSE, shrinkage==FALSE, boot==0))) {
    Result=sveig_fastsvd(data_matrix)
  } else if (boot>0) {

    if (rarefy>1) { rarefy = rarefy} else { rarefy = nrow(data_matrix) }

    if (shrinkage==FALSE) {
      sveig_boot=do.call("rbind", lapply(seq(boot), function(x) {
        indices_boot=sample(seq(nrow(data_matrix)), rarefy, replace = TRUE)
        return(sveig_fastsvd(data_matrix[indices_boot,]))
          }))
    }

    if (shrinkage==TRUE) {
      sveig_boot=do.call("rbind", lapply(seq(boot), function(x) {
        indices_boot=sample(seq(nrow(data_matrix)), rarefy, replace = TRUE)
        cov_mat=nlshrink::linshrink_cov(data_matrix[indices_boot,])
        return(sveig_covmat(cov_mat))
          }))
      }

    boot_ci_sveig=c(sort(sveig_boot[,1])[c(round(boot*0.025), round(boot*0.975))],
                    sort(sveig_boot[,2])[c(round(boot*0.025), round(boot*0.975))]
      )

    Result=list(mean_boot_sveig=mean(sveig_boot[,1]),
                median_dimensions_used=median(sveig_boot[,2]),
                boot_95ci_sveig=boot_ci_sveig,
                boot_sveig=sveig_boot)

  } else if (boot==0) {
    if (rarefy==FALSE & shrinkage==TRUE) {
      cov_mat=nlshrink::linshrink_cov(data_matrix)
      Result=sveig_covmat(cov_mat)
    } else if (rarefy>1) {
      warning("Only a SINGLE rarefaction will be performed.
 This estimate alone is unreliable.
 Please, repeat the analysis multiple times and compute the mean to obtain sensible estimates")
      if (shrinkage==FALSE) {

      indices_rar=sample(seq(nrow(data_matrix)), rarefy, replace = FALSE)
      Result=sveig_fastsvd(data_matrix[indices_rar,])

      } else if (shrinkage==TRUE) {
          indices_rar=sample(seq(nrow(data_matrix)), rarefy, replace = FALSE)
          cov_mat=nlshrink::linshrink_cov(data_matrix[indices_rar,])
        Result=sveig_covmat(cov_mat)
        }

      }
  }
return(Result)
}


# Function to use fast.svd from corpcor to compute
# positive eigenvalues from a data matrix
# and then compute relative variance of eigenvalues
sveig_fastsvd=function(data_matrix) {
  sample_size=nrow(data_matrix)
  eigenvalues=(corpcor::fast.svd(
  scale(data_matrix, center = TRUE, scale=FALSE))$d^2)/(sample_size-1)
  eigenvalues=eigenvalues/sum(eigenvalues)
  veig=var_pop(eigenvalues)
  dimensions=length(eigenvalues)
  max_theoric_veig=(dimensions-1)/(dimensions^2)
  scal_var_eig=veig/max_theoric_veig
  res=c(scal_var_eig, dimensions)
  names(res)=c("Scaled_variance_of_eigenvalues", "dimensions")
return(res)
}

# Function to compute an eigendecomposition of a covariance matrix
# and then compute relative variance of eigenvalues
sveig_covmat=function(covmat) {
  eigenvalues=eigen(covmat)$values
  eigenvalues=eigenvalues[eigenvalues>0]
  eigenvalues=eigenvalues/sum(eigenvalues)
  dimensions=length(eigenvalues)
  veig=var_pop(eigenvalues)
  max_theoric_veig=(dimensions-1)/(dimensions^2)
  scal_var_eig=veig/max_theoric_veig
  res=c(scal_var_eig, dimensions)
  names(res)=c("Scaled_variance_of_eigenvalues", "dimensions")
  return(res)
return(res)
}

# Function to compute population-level variance
# (by dividing the sum of squares for the number of observations
# rather than the number of observations minus one)
var_pop=function(x) {
  sqdev=(x-mean(x))^2
  sum(sqdev)/length(x)
}
