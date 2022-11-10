#' Bootstrapped distance between two arrays
#'
#' Computes bootstrapped estimates of the mean distance
#' between two groups and their confidence intervals.
#'
#' This may be useful to compare whether the differences between two groups
#' are larger or smaller than differences between two other groups.
#'
#' For instance, if we wanted to quantify shape sexual dimorphism in two populations,
#' we could run this analysis separately for the two populations and then check the
#' confidence intervals. If the two confidence intervals are disjunct there is evidence
#' for the two populations having different levels of sexual dimorphism.
#'
#' The computation performs bootstrap by resampling with replacement
#' within each of the two groups and at each round computing the Euclidean distance
#' between the two groups.
#' It is also possible to resample at a different sample size than the one in the data
#' using the attributes nA and nB.
#'
#' Notice that the confidence interval is expressed on a scale between 0 and 1 and not
#' as a percentage (e.g., 0.95 means 95% confidence interval)
#'
#'
#' @param A,B Matrices or data frames containing data
#' (observations in rows, variables in columns).
#' @param boot number of bootstrap resamples
#' @param ci width of the confidence interval
#' @param nA,nB sample sizes for each bootstrapped group (defaults to original sample size)
#'
#' @return The function outputs a named vector with the mean, median, upper and lower confidence interval bounds
#' obtained from the bootstrapped samples
#'
#'
#' @import stats
#' @export

dist_mean_boot=function(A, B, boot=1000, ci=0.95, nA=nrow(A), nB=nrow(B)){
  if (any(ci<=0, ci>=1)){stop("The width of the confidence interval should be higher than 0 and lower than 1")}
  if (any(nrow(A)<=1, nrow(B)<=1)){stop("Both matrices should have more than one row and more than one column")}
  if (any(ncol(A)<=1, ncol(B)<=1)){stop("Both matrices should have more than one row and more than one column")}
  if (ncol(A)!=ncol(B)){stop("The number of columns of the two matrices is not the same, cannot compute distances")}

  ssA=nrow(A) ; ssB=nrow(B)

  boot_idx=lapply(seq(boot), function(i)
                  list(Aboot=sample(seq(ssA), nA, replace = TRUE),
                       Bboot=sample(seq(ssB), nB, replace = TRUE)))
  boot_dist=unlist(lapply(seq(boot), function(i){
    dist(rbind(colMeans(A[boot_idx[[i]]$Aboot,]),
         colMeans(B[boot_idx[[i]]$Bboot,])))
    }))
  ci_bnd_idx=c(round(boot*((1-ci)/2)), round(boot*(((1-ci)/2)+ci)))
  boot_ci=sort(boot_dist)[ci_bnd_idx]
  results=c(mean(boot_dist), median(boot_dist), boot_ci)
  names(results)=c("Mean distance", "Median distance", "Lower CI extreme", "Upper CI extreme")
return(results)
}
