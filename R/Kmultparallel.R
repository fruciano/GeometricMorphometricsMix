#' Parallel implementation of Adams' Kmult
#'
#' Parallel implementation of Kmult, a measure of phylogenetic signal which
#' is a multivariate equivalent of Blomberg's K
#'
#'
#'  This is an updated and improved version of the function included in Fruciano et al. 2017.
#'  It performs the computation of Adams' Kmult (Adams 2014) in parallel
#'  with the aim of facilitating computation on a distribution of trees rather than a single tree
#'  This is a very simple parallelization of the code ('embarassingly parallel') which
#'  uses internally Adams' published code (to maintain consistency over time).
#'  If one wanted to perform a computation of Kmult on a single tree, he/she would be
#'  advised to use the version implemented in the package geomorph.
#'
#'
#' @section Notice:
#' Under Windows, this implementation needs the package 'parallelsugar'
#'if one wants to run it with multiple cores in parallel.
#' The package 'parallelsugar' can be installed from Github with the command
#' devtools::install_github("nathanvan/parallelsugar")
#' This, in turn, requires the package devtools to be installed.
#' Notice also that the reason why the value of iter is set to 0 by default is that
#' using permutations slows down sensibly the computation and it is of doubtful utility
#' in most cases (as one would get a distribution of p-values not independent from each other).
#' Unless one has a specific use in mind, the suggestion is, then, to keep iter at the default value
#'
#'
#' @section Citation:
#' If you use this function please kindly cite both
#' Fruciano et al. 2017 (because you're using this parallelized function)
#' Adams 2014 (because the function computes Adams' Kmult)
#'
#' @param data data frame with continuous (multivariate) phenotypes,
#' including geometric morphometric data
#' @param trees a multi-phylo object containing a collection of trees
#' @param ncores number of cores to use (default 1; i.e., no parallelization)
#' @param burninpercent percentage of trees in trees to discard as burn-in
#' (by default no tree is discarded)
#' @param iter number of permutations to be used in the permutation test
#' (this should normally be left at the default value of 0)
#'
#' @return The function outputs a matrix with as many rows as input trees and columns:
#'  \describe{
#'   \item{Phylogenetic signal (K mult)}{Value of Kmult for each tree}
#'   \item{p value}{p value for the significance of the test (normally not used)}
#' }
#'
#' @references Adams DC. 2014. A Generalized K Statistic for Estimating Phylogenetic Signal from Shape and Other High-Dimensional Multivariate Data. Systematic Biology 63:685-697.
#' @references Fruciano C, Celik MA, Butler K, Dooley T, Weisbecker V, Phillips MJ. 2017. Sharing is caring? Measurement error and the issues arising from combining 3D morphometric datasets. Ecology and Evolution 7:7034-7046.
#'
#' @export
Kmultparallel = function(data, trees, ncores = 1, burninpercent = 0, iter = 0) {
    trees = trees[round((length(trees)/100) * burninpercent):length(trees)]
    # Removes the percentage of trees to be considered as burn-in (rounded to the next integer)

    droplist = setdiff(trees[[1]]$tip, row.names(data))
    # Use the names of the tips of the first tree and the row names of the data to determine which tips should be pruned
    # (i.e., it assumes that all trees share the same tips and these are in the same order)


    pruning = function(tree) {
        prunedtree = ape::drop.tip(tree, droplist)
        return(prunedtree)
    }
	
	# If a single core is used, use lapply,
	# otherwise use mclapply (possibly with the parallelsugar implementation)
	if (ncores==1) {
					mclapply=lapply
					} else {	
    if (Sys.info()["sysname"] == "Windows") {	
        if (!requireNamespace("parallelsugar", quietly = T)) {
		stop("Package \"parallelsugar\" needed for this function to work under Windows.\names
		Please install it from Github devtools::install_github(\"nathanvan/parallelsugar\") .",
		call. = FALSE)
        }
     }

    mclapply <- switch(Sys.info()[["sysname"]], Windows = {
        parallelsugar::mclapply_socket
    }, Linux = {
        parallel::mclapply
    }, Darwin = {
        parallel::mclapply
    })
	 
	} 
	prunedtrees = mclapply(trees, pruning, mc.cores = ncores)
    class(prunedtrees) = "multiPhylo"
    # Use mclapply for speeding up the pruning of the trees Note that under Windows, it will install the package parallelsugar
    # (using devtools, which must be installed) as otherwise mclapply does not work with multiple cores

    datareordered = mclapply(prunedtrees, function(x) data[x$tip, ])
    # Reorder data with the same order of the tips

    ntrees = seq_len(length(prunedtrees))
    kmulttrees = mclapply(ntrees, mc.cores = ncores, function(x) Test_Kmult(datareordered[[x]], prunedtrees[[x]], iter = iter))

    kmulttreesmatrix = t(matrix(simplify2array(unlist(kmulttrees, recursive = F)), 2))
    rownames(kmulttreesmatrix) = names(prunedtrees)
    colnames(kmulttreesmatrix) = c("Phylogenetic signal (K mult)", "p value")
    # Produce a matrix with the results, which gets returned to the user
    if (iter == 0) {
        return(cbind(kmulttreesmatrix[, 1]))
    } else {
        return(kmulttreesmatrix)
    }
}



### Function to compute Kmult, as published by Adams (not exported)
Test_Kmult <- function(x, phy, iter = 999) {
    Kmult <- function(x, phy) {
        x <- as.matrix(x)
        N <- length(phy$tip.label)
        ones <- array(1, N)
        C <- ape::vcv.phylo(phy)
        C <- C[row.names(x), row.names(x)]
        a.obs <- colSums(solve(C)) %*% x/sum(solve(C))  #evol.vcv code
        distmat <- as.matrix(dist(rbind(as.matrix(x), a.obs)))
        MSEobs.d <- sum(distmat[(1:N), (N + 1)]^2)  #sum distances root vs. tips
        eigC <- eigen(C)
        D.mat <- solve(eigC$vectors %*% diag(sqrt(eigC$values)) %*% t(eigC$vectors))
        dist.adj <- as.matrix(dist(rbind((D.mat %*% (x - (ones %*% a.obs))), 0)))
        MSE.d <- sum(dist.adj[(1:N), (N + 1)]^2)  #sum distances for transformed data)
        K.denom <- (sum(diag(C)) - N * solve(t(ones) %*% solve(C) %*% ones))/(N - 1)
        K.stat <- (MSEobs.d/MSE.d)/K.denom
        return(K.stat)
    }
    K.obs <- Kmult(x, phy)

    P.val <- 1
    K.val <- rep(0, iter)
    for (i in 1:iter) {
        x.r <- as.matrix(x[sample(nrow(x)), ])
        rownames(x.r) <- rownames(x)
        K.rand <- Kmult(x.r, phy)
        P.val <- ifelse(K.rand >= K.obs, P.val + 1, P.val)
        K.val[i] <- K.rand
    }
    P.val <- P.val/(iter + 1)
    K.val[iter + 1] = K.obs
    hist(K.val, 30, freq = TRUE, col = "gray", xlab = "Phylogenetic Signal")
    arrows(K.obs, 50, K.obs, 5, length = 0.1, lwd = 2)
    return(list(phy.signal = K.obs, pvalue = P.val))
}
