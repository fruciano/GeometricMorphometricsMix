### Please, if you use this function, kindly cite:
### Fruciano et al. 2017 - Ecology and Evolution (because you're using this parallelized function)
### Adams 2014 - Systematic Biology (because the function computes Adams' Kmult)

### This is an updated version of the function included in Fruciano et al. 2017
### This function performs the computation of Adams' Kmult (Adams 2014 - Systematic Biology) in parallel
### with the aim of facilitating computation on a distribution of trees rather than a single tree
### This is a very simple parallelization of the code ("embarassingly parallel")
### The function assumes that the code published in Adams' paper (Appendix 2) has been saved in the file "Test.Kmult.R"
### and that this file is located in the working directory (as it is "sourced"); this choice has been made to improve usability of code over time
### Differently from the version in the supplementary of Fruciano et al. 2017,
### this version DOES NOT assume that the trees in the distribution have all the same tips and these are ordered consistently
### Alternative implementations using the package geomorph or performing the computation in a serial fashion can be requested to Carmelo Fruciano
### Under Windows, this implementation installs (if not yet installed) the package "parallelsugar" from Github
### To this aim, devtools should be installed prior to using the function

### Usage: Kmultparallel(data,trees,ncores=1,burninpercent=0,iter=0)
### data is a data frame with continuous (multivariate) phenotypes, including geometric morphometric data
### trees is a multi-phylo object containing the trees
### ncores is the number of cores to use (default 1; i.e., no parallelization)
### burninpercent is the percentage of trees in trees to discard as burn-in (by default no tree is discarded)
### iter is the number of permutations to be used in the analysis of each Kmult value (this is set by default to 0
### as it is normally not useful to obtain highly dependent p-values from a distribution of trees)

### The function outputs a matrix with the value of Kmult and the p-value for each of the trees in trees

Kmultparallel=function(data,trees,ncores=1,burninpercent=0,iter=0) {
library(ape)
library(parallel)
trees=trees[round((length(trees)/100)*burninpercent):length(trees)]
# Removes the percentage of trees to be considered as burn-in (rounded to the next integer)

droplist=setdiff(trees[[1]]$tip, row.names(data))
# Use the names of the tips of the first tree and the row names of
# the data to determine which tips should be pruned
# (i.e., it assumes that all trees share the same tips and these are in the same order)


pruning=function(tree) {
prunedtree=drop.tip(tree,droplist)
return(prunedtree)
}


if(Sys.info()['sysname']=="Windows") {

if (!require("parallelsugar", character.only=T, quietly=T)) {
library(devtools)
install_github('nathanvan/parallelsugar')
}
library(parallelsugar)
prunedtrees=mclapply_socket(trees, pruning, mc.cores = ncores)

 } else { 
prunedtrees=mclapply(trees, pruning, mc.cores = ncores)
 }
class(prunedtrees)="multiPhylo"
# Use mclapply for speeding up the pruning of the trees
# Note that under Windows, it will install the package parallelsugar
# (using devtools, which must be installed) as otherwise mclapply
# does not work with multiple cores

mclapply <- switch( Sys.info()[['sysname']],
                     Windows = {mclapply_socket},
                     Linux   = {mclapply},
                     Darwin  = {mclapply})


datareordered=mclapply(prunedtrees, function(x) data[x$tip,])
# Reorder data with the same order of the tips

source('Test.Kmult.R')

ntrees=1:length(prunedtrees)
kmulttrees=mclapply(ntrees, mc.cores=ncores, function(x) Test.Kmult(datareordered[[x]],prunedtrees[[x]],iter=iter))

kmulttreesmatrix=t(matrix(simplify2array(unlist(kmulttrees, recursive=F)), 2))
rownames(kmulttreesmatrix)=names(prunedtrees)
colnames(kmulttreesmatrix)=c("Phylogenetic signal (K mult)", "p value")
# Produce a matrix with the results, which gets returned to the user
return(kmulttreesmatrix)
}

