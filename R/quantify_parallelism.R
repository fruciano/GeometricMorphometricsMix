#' Quantify and test for parallelism in multivariate vectors
#'
#' @description
#' Computes various metrics to quantify parallelism among multivariate vectors,
#' including vector correlations, eigenvalue analysis, and statistical tests.
#'
#' @details
#' The function largely follows the approaches outlined in De Lisle & Bolnick (2020)
#' and Watanabe (2022), whereas the Cramer's V-type scaling
#' for the Rayleigh test statistic is presented in Fruciano, Franchini et al. in preparation.
#' Specifically, the Schott test (Schott 2005) and the Rayleigh test (Mardia & Jupp 1999)
#' are discussed estensively in Watanabe (2022).
#'
#' @param X A matrix where each row represents a multivariate vector
#' @param test Logical. If TRUE, performs statistical tests for vector concentration (default: TRUE)
#' @param ndim Number of dimensions in the vector space (default: number of columns in X)
#'
#' @return A list containing:
#' - `C`: Correlation matrix
#' - `eigen_C`: Eigendecomposition of the correlation matrix
#' - `A`: Transformation matrix
#' - `prop_v`: Proportion of variance explained by first eigenvalue
#' - `var_eigenval_C`: Variance of eigenvalues
#' - `Schott_test`: Schott's test results for vector concentration
#' - `Rayleigh_test`: Rayleigh test results for vector concentration
#'
#' @references De Lisle, S. P., Bolnick, D.I. 2020. A multivariate view of parallel evolution. Evolution 74(7): 1466-1481.
#' @references Watanabe, J. 2022. Detecting (non)parallel evolution in multidimensional spaces: angles, correlations and eigenanalysis. Biology Letters 18(2):20210638.
#' @references Fruciano, Franchini et al. in preparation. Limited genetic parallelism underlies repeated sympatric divergence in Nicaraguan cichlid fish.
#' @references Mardia, K. V. and P. E. Jupp (1999). Directional Statistics.
#' @references Schott, J. R. 2005. Testing for complete independence in high dimensions. Biometrika 92(4): 951-956.
#'
#' @examples
#'
#' # In this example we will create multivariate data for four random groups
#' # (we could imagine this is morphometric data for four populations)
#' # and we will test whether the difference between groups A and B is parallel
#' # to the difference between groups C and D.
#'
#' set.seed(123)
#' library(MASS)
#' X=as.data.frame(mvrnorm(100, mu=rep(0, 50), Sigma=diag(50)))
#' # Create a matrix of 100 observations from a multivariate normal distribution
#' split_fac=as.factor(rep(c("A", "B", "C", "D"), each=25))
#' # Create a factor with 4 levels representing 4 hypothetical groups
#' X_split=split(X, f=split_fac)
#' # split the dataset based on the factor
#' X_means=do.call("rbind", lapply(X_split, colMeans))
#' # Compute vectors of average for each group
#' # (notice that, as we obtained random data and we are splitting into
#' # 4 random groups, we do not expect parallelism between the vectors)
#' #
#' # Now let us compute vectors of differences between the average of groups in pairs
#' X_diffs=rbind((X_means[1,]-X_means[2,]), (X_means[3,]-X_means[4,]))
#' # Notice that each row is a vector of differences
#'
#' QP=quantify_parallelism(X_diffs)
#' # Run the function to quantify parallelism
#' # Checking the results, it will become clear how the vectors are not parallel
#' # This can be seen, for example, by looking at the p values for the Schott and Rayleigh tests
#'
#' @export
quantify_parallelism=function(X, test=TRUE, ndim=ncol(X)){
  # Notation follows the paper as much as possible
    X=t(apply(X,1, unit_vector_scale))
  # Scale vectors to unit size in case they are not scaled

  C=X%*%t(X)
  eigen_C=eigen(C)
  # Compute the eigendecomposition of C following De Lisle & Bolnick 2020 - Evolution

  Results=list(C=C,
               eigen_C=list(Q=eigen_C$vector,
                            V=eigen_C$values),
               A=t(X)%*%solve(t(eigen_C$vector)),
               prop_v=eigen_C$values[1]/sum(eigen_C$values),
               var_eigenval_C=var(eigen_C$values))
  if (test==TRUE){
    # Commenting out this part because, as of 14/11/2024 the function
    # Directional::rvmf() called by function sim_sumsqvecor()
    # stopped working with attribute k=0 (uniform spherical distribution of vectors)
    # For this reason, all tests will be done without the simulation version
    # (earlier tests were done with the simulation as it was working)
    # if (n_sim_test>0){
    #  simulat_veccor=sim_sumsqvecor(n=nrow(X), p=ncol(X), nsim=n_sim_test)
    #  Results$Schott_test=schott_test(X, C,
    #                           expct_sumsqrcor=simulat_veccor[1],
    #                           var_sumsqrcor=simulat_veccor[2], p=ndim)
    #  Results$Rayleigh_test=Rayleigh(X, p=ndim, nsim=n_sim_test)
    # } else {
    Results$Schott_test=schott_test(X, C, p=ndim)
    Results$Rayleigh_test=Rayleigh(X, p=ndim)
    # }
  }
return(Results)
}


unit_vector_scale=function(vec) {vec / sqrt(sum(vec^2))}

schott_test=function(X, C, expct_sumsqrcor=NULL, var_sumsqrcor=NULL, p){
  # Schott's (2005) for independent directions
  # as per Watanabe 2022 - Biology Letters
  # (significant values indicate concentration of vectors)

  n=nrow(X)
  # Number of vectors of differences, and number of variables

  veccor=cor(t(X))
  sumsqvecor=sum(veccor[lower.tri(veccor, diag=FALSE)]^2)
  # Sum of squared vector correlations

  if (all(is.null(expct_sumsqrcor), is.null(var_sumsqrcor))){
    exp_sumsqvecor=(n*(n-1))/(2*p)
    exp_var_smsqvecor=(n*(n-1)*(p-1))/(p^2*(p+2))
    # Formulas 2.16 in Watanabe 2022

    Z_score = (sumsqvecor - exp_sumsqvecor)/sqrt(exp_var_smsqvecor)
  } else {
    Z_score = (sumsqvecor - expct_sumsqrcor)/sqrt(var_sumsqrcor)
  }

    p_value=pnorm(Z_score, lower.tail = FALSE)

return(data.frame(Z_score=Z_score, p_value=p_value, row.names = NULL))
}

# Commenting this out as there are problems with Directional::rvmf() as of 14/11/2024
# (see above). When these will be fixed, remember to import rvmf from the package
# instead of using the call to Directional::rvmf()
# sim_sumsqvecor=function(n, p, nsim=1000){
#   vecs=lapply(seq(nsim), function(i) {
#     Directional::rvmf(n=n, mu=unit_vector_scale(rnorm(p)), k=0)
#     })
#   vecs=lapply(vecs, function(X) t(apply(X, 1, unit_vector_scale)))
#   # Random unit length vectors
#
#   sumsqrvecor=unlist(lapply(vecs, function(Y){
#     veccor=cor(t(Y))
#     sumsqvecor=sum(veccor[lower.tri(veccor, diag=FALSE)]^2)
#   }))
#   # Compute sum of squared vector correlations for the random vectors
#
#   res=c(mean=mean(sumsqrvecor), variance=var(sumsqrvecor))
# return(res)
# }
#
rayleigh_stat=function(X, p, n){
  zbar=colMeans(X)
  S=n*p*crossprod(zbar)
  # Formula B.1 in Watanabe 2022 - Biology Letters
  Sstar=((1-(1/(2*n)))*S)+((1/(2*n*(p+2)))*(S^2))
  # Formula B.2 in Watanabe 2022 - Biology Letters
  return(c(S=S, Sstar=Sstar))
}


rayleigh_pval=function(S, Sstar, p, n, nsim=0){
  if (nsim>0){
    sim_vecs=lapply(seq(nsim), function(i){
      Directional::rvmf(n=n, mu=rep(0, p), k=0)
    })
    S_dist=do.call("rbind", lapply(sim_vecs, function(Y){
      rayleigh_stat(Y, p, n)
    }))
    pval_S=(sum(S_dist[,1]>=S)+1)/(nsim+1)
    pval_Sstar=(sum(S_dist[,2]>=Sstar)+1)/(nsim+1)
  }else{
    pval_S=pchisq(S, p, lower.tail = FALSE)
    pval_Sstar=pchisq(Sstar, p, lower.tail = FALSE)
  }
return(c(pval_S=pval_S, pval_Sstar=pval_Sstar))
}

Vscale=function(x, n, df) {sqrt(x/(n*df))}
# To obtain a Cramer’s V type scaling to that the statistic of the test is scaled
# by the number of variables

Rayleigh=function(X, p, nsim=0){
  n=nrow(X)
  rstat=rayleigh_stat(X, p, n)
  rstat_scaled=unlist(lapply(rstat, function(x) Vscale(x, df=p, n=n)))
  names(rstat_scaled)=c("S_Vscaled", "Sstar_Vscaled")
  pvalues=rayleigh_pval(S=rstat[1], Sstar=rstat[2], p, n, nsim)
  return(c(rstat, rstat_scaled, pvalues))
}
