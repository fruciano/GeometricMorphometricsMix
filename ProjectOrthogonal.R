### This function projects data to the subspace orthogonal to a multivariate column vector
### Both data and vector should be in matrix form

### This is the approach used in Burnaby's (1966 - Biometrics) method for size correction
### and in Valentin et al. (2008 - Journal of Fish Biology) to remove the effect of dorso-ventral arching in fish
### (see also Fruciano 2016 - Development Genes and Evolution for a discussion and other examples of use)

### Usage: ProjectOrthogonal(Data,vector)
### Data is the original data matrix
### vector is the vector that should be swept from the data (provided as a matrix with one column)

### The notation here follows as much as possible Rohlf & Bookstein (1987 - Systematic Zoology) for ease of reference

ProjectOrthogonal=function(Data,vector) {
Ip=diag(nrow(vector))
F=vector
L=Ip-F%*%(solve(t(F)%*%F))%*%t(F)
ProjectedData=Data%*%L
return(ProjectedData)
}

