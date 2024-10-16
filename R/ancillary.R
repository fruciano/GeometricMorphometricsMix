### Ancillary functions
### (not exported)


# Functions to repeat rows or columns
#' @noRd
rep.row=function(x,n){
   matrix(rep(x,each=n),nrow=n)
}
#' @noRd
rep.col=function(x,n){
   matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}


# Functions to convert radians to degrees
# and degrees to radians
rad2deg = function(rad) {(rad * 180) / (pi)}
deg2rad = function(deg) {(deg * pi) / (180)}


# Function to compute "multivariate variance"
# (sum of univariate variances)
muvar=function(X) {
  sum(apply(X, 2, var))
}

# Function to compute "Proper variance" as described in
# Claramunt 2010 - Evolution
# It uses a linear shrinkage estimator of the covariance matrix
claramunt_proper_variance=function(X) {
  covX=nlshrink::linshrink_cov(as.matrix(X))
  eigenv=eigen(covX)$values
  proper_variance=(sum(sqrt(eigenv)))^2
return(proper_variance)
}


# Function to compute mean pairwise Euclidean distances
meanpairwiseEuclideanD=function(X) {
  mean(dist(X,method = "euclidean"))
}

# Function to compute the Euclidean distance
# between two observations
pdistance=function(X1, X2) {
  dist(rbind(X1, X2),method = "euclidean")
}

