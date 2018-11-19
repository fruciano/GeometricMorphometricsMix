### Ancillary functions
### (not exported)


# Functions to repeat rows or columns
rep.row=function(x,n){
   matrix(rep(x,each=n),nrow=n)
}
rep.col=function(x,n){
   matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}


# Functions to convert radians to degrees
# and degrees to radians
rad2deg = function(rad) {(rad * 180) / (pi)}
deg2rad = function(deg) {(deg * pi) / (180)}


# Function to compute "multivariate variance"
# (trace of the covariance matrix)
muvar=function(X) {
  sum(diag(cov(X)))
}

# Function to compute mean pairwise Euclidean distances
meanpairwiseEuclideanD=function(X) {
  mean(dist(X,method = "euclidean"))
}

