% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Kmultparallel.R
\name{Kmultparallel}
\alias{Kmultparallel}
\title{Parallel implementation of Adams' Kmult}
\usage{
Kmultparallel(data, trees, ncores = 1, burninpercent = 0, iter = 0)
}
\arguments{
\item{data}{data frame with continuous (multivariate) phenotypes,
including geometric morphometric data}

\item{trees}{a multi-phylo object containing a collection of trees}

\item{ncores}{number of cores to use (default 1; i.e., no parallelization)}

\item{burninpercent}{percentage of trees in trees to discard as burn-in
(by default no tree is discarded)}

\item{iter}{number of permutations to be used in the permutation test
(this should normally be left at the default value of 0)}
}
\value{
The function outputs a matrix with as many rows as input trees and columns:
 \describe{
  \item{Phylogenetic signal (K mult)}{Value of Kmult for each tree}
  \item{p value}{p value for the significance of the test (normally not used)}
}
}
\description{
Parallel implementation of Kmult, a measure of phylogenetic signal which
is a multivariate equivalent of Blomberg's K
}
\details{
This is an updated and improved version of the function included in Fruciano et al. 2017.
 It performs the computation of Adams' Kmult (Adams 2014) in parallel
 with the aim of facilitating computation on a distribution of trees rather than a single tree
 This is a very simple parallelization of the code ('embarassingly parallel') which
 uses internally Adams' published code (to maintain consistency over time).
 If one wanted to perform a computation of Kmult on a single tree, he/she would be
 advised to use the version implemented in the package geomorph.
}
\section{Notice}{

Under Windows, this implementation needs the package 'parallelsugar'
if one wants to run it with multiple cores in parallel.
The package 'parallelsugar' can be installed from Github with the command
devtools::install_github("nathanvan/parallelsugar")
This, in turn, requires the package devtools to be installed.
Notice also that the reason why the value of iter is set to 0 by default is that
using permutations slows down sensibly the computation and it is of doubtful utility
in most cases (as one would get a distribution of p-values not independent from each other).
Unless one has a specific use in mind, the suggestion is, then, to keep iter at the default value
}

\section{Citation}{

If you use this function please kindly cite both
Fruciano et al. 2017 (because you're using this parallelized function)
Adams 2014 (because the function computes Adams' Kmult)
}

\references{
Adams DC. 2014. A Generalized K Statistic for Estimating Phylogenetic Signal from Shape and Other High-Dimensional Multivariate Data. Systematic Biology 63:685-697.

Fruciano C, Celik MA, Butler K, Dooley T, Weisbecker V, Phillips MJ. 2017. Sharing is caring? Measurement error and the issues arising from combining 3D morphometric datasets. Ecology and Evolution 7:7034-7046.
}
