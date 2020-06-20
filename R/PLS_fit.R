#' Partial least squares (PLS) analysis
#'
#' Performs a two-block PLS analysis, optionally allowing
#' for tests of significance using permutations
#'
#' This function performs a PLS analysis (sensu Rohlf & Corti 2000).
#' Given two blocks of variables (shape or other variables) scored on the same observations (specimens),
#' this analysis finds a series of pairs of axis accounting for maximal covariance between the two blocks.
#' If tests of significance with permutations are selected, three different significance tests are performed:
#'  \itemize{
#'   \item{Global significance:}{ tested using Escoufier RV}
#'   \item{Axis-specific significance based on singular value:}{ this is the same test described in Rohlf & Corti 2000}
#'   \item{Axis-specific significance based on correlation of PLS scores:}{ this is a commonly used test which uses as statistic, for each pair of PLS (singular) axes, the correlation of the scores of the first block with the scores of the second block}
#' }
#' The object of class pls_fit returned by the function has print() and summary() methods associated to it.
#' This means that using these generic functions on an object created by this function (see examples), it is possible to obtain information on the results.
#' In particular, print() returns a more basic set of results on the global association, whereas summary() returns (only if permutation tests are used) results for each pair of singular axes.
#'
#'
#' @section Notice:
#' \itemize{
#' \item{The function does NOT perform GPA when applied to separate configurations of points.}
#' \item{When using the Escoufier RV, notice that the value reported is the observed value without rarefaction.
#' For a description of the problem, please see Fruciano et al 2013. To obtain rarefied estimates of Escoufier RV and their confidence interval, use the function RVrarefied.}
#' \item{In the permutation test, rows of Y are permuted, so using the block with fewer variables as Y may speed up computations and substantially reduce memory usage}
#' \item{When using the print() and summary() on the pls_fit objects obtained with this function, some of the values are rounded for ease of interpretation. The non-rounded values can be obtained accessing individual elements of the object (see examples).}
#'}
#'
#' @section Citation:
#' If you use this function to perform the PLS analysis and test for significance,
#' cite Rohlf & Corti 2000 (or earlier references outside of geometric morphometrics).
#' If you report the test of significance based on the Escoufier RV coefficient, also cite Escoufier 1975.
#' If you also use the predict() method to obtain estimates of the shape (or other variable)
#' predicted by each pair of axis, please cite Fruciano et al. under review
#' (contact Carmelo Fruciano for the most recent reference)
#'
#' @param X,Y Matrices or data frames
#' containing each block of variables
#' (observations in rows, variables in columns).
#' @param perm number of permutations to use for hypothesis testing
#'
#' @seealso \code{\link{RVrarefied}}
#' @return The function outputs an object of class "pls_fit" and "list" with the following elements:
#'  \describe{
#'   \item{XScores}{Scores along each singular (PLS) axis for the first block of variables (X)}
#'   \item{YScores}{Scores along each singular (PLS) axis for the first block of variables (Y)}
#'   \item{U}{Left singular axes}
#'   \item{V}{Right singular axes}
#'   \item{D}{Singular values}
#'   \item{percentage_squared_covariance}{Percented squared covariance accounted by each pair of axes}
#'   \item{global_significance_RV}{(only if perm>0) Observed value of Escoufier RV coefficient and p value obtained from the permutation test}
#'   \item{singular_axis_significance}{(only if perm>0) For each pair of singular (PLS) axis, the singular value, the correlation between scores, their significance level based on permutation, and the proportion of squared covariance accounted are reported}
#'   \item{OriginalData}{Data used in the analysis}
#'   \item{x_center}{Values used to center data in the X block}
#'   \item{y_center}{Values used to center data in the Y block}
#'
#' }
#'
#'
#'
#' @references Escoufier Y. 1973. Le Traitement des Variables Vectorielles. Biometrics 29:751-760.
#' @references Fruciano C, Colangelo P, Castiglia R, Franchini P. under review. Does divergence from normal patterns of integration increase as chromosomal fusions increase in number? A test on a house mouse hybrid zone.
#' @references Rohlf FJ, Corti M. 2000. Use of Two-Block Partial Least-Squares to Study Covariation in Shape. Systematic Biology 49:740-753.
#' @examples
#'
#' ##############################
#' ### Example 1: random data ###
#' ##############################
#'
#' library(MASS)
#' set.seed(123)
#' A=mvrnorm(100,mu=rep(0,20), Sigma=diag(20))
#' B=mvrnorm(100,mu=rep(0,10), Sigma=diag(10))
#' # Create two blocks of, respectively, 20 and 10 variables
#' # for 100 observations.
#' # This simulates two different blocks of data (shape or otherwise) measured on the same individuals.
#' # Note that, as we are simulating them independently,
#' # we don't expect substantial covariation between blocks
#'
#' PLS_AB=pls(A, B, perm=99)
#' # Perform PLS analysis and use 99 permutations for testing
#' # (notice that in a real analysis, normally one uses more permutations)
#' print(PLS_AB)
#' # As expected, we do not find significant covariation between the two blocks
#'
#' summary(PLS_AB)
#' # The same happens when we look at the results for each of the axes
#'
#' # Notice that both for print() and summary() some values are rounded for ease of visualization
#' # However, the correct values can be always obtained from the object created by the function
#' # e.g.,
#' PLS_AB$singular_axis_significance
#'
#'
#' ######################################
#' ### Example 2: using the classical ###
#' ### iris data set as a toy example ###
#' ######################################
#'
#' data(iris)
#'# Import the iris dataset
#' set.seed(123)
#'
#'versicolor_data=iris[iris$Species=="versicolor",]
#'# Select only the specimens belonging to the species Iris versicolor
#'versicolor_sepal=versicolor_data[,grep("Sepal", colnames(versicolor_data))]
#'versicolor_petal=versicolor_data[,grep("Petal", colnames(versicolor_data))]
#'# Separate sepal and petal data
#'PLS_sepal_petal=pls(versicolor_sepal, versicolor_petal, perm=99)
#'# Perform PLS with permutation test
#'# (again, chosen few permutations)
#'
#'print(PLS_sepal_petal)
#'summary(PLS_sepal_petal)
#'# Global results and results for each axis (suggesting significant association)
#'
#'
#' @export
pls=function(X, Y, perm=999) {
  if ((!(class(X) %in% c("matrix", "data.frame"))) |  (!(class(Y) %in% c("matrix", "data.frame")))) {
    stop("X and Y should be matrices or data.frames (observations in rows, variables in columns)")
  }
  if (nrow(X)!=nrow(Y)) {
    stop("X and Y should have the same number of rows (different variables/shape scored on the same individuals)")
  }
  # Checks of class and consistent number of rows

  if (perm>0) {
    Results=pls_perm(X, Y, perm=perm)
  } else {
    Results=pls_base(X, Y)
  }
return(Results)
}
# Wrapper for using pls_base or pls_perm, depending on the presence of permutations


#' @export
print.pls_fit = function (x, ...) {
  cat("Two-block PLS analysis")
  cat("\n")
  cat("======================")
  cat("\n")
  if ("global_significance_RV" %in% names(x)) {
   cat("Global significance based on permuting the Escoufier RV coefficient")
    cat("\n")
    cat(paste0("Observed RV=", round(x$global_significance_RV["observed_RV"], 2)))
    cat("\n")
    cat(paste0("p value= ", x$global_significance_RV["RV_p_value"]))
    cat("\n")
    cat("Note that the observed RV is not rarefied")
    cat("\n")
    cat("\n")
    cat("For axis-specific results, please use summary()")
    cat("\n")
  } else {
    Res=data.frame(Axis=seq(length(x$D)),
                   percentage_squared_covariance=x$percentage_squared_covariance,
                   singular_value=x$D,
                   correlation_between_scores=x$CorrXScoresYScores)
    print(Res)
  }

  cat("======================")

}

# print method for S3 class pls_fit
#
#' @export
summary.pls_fit = function(object, ...) {
  if ("global_significance_RV" %in% names(object)) {
  cat("Two-block PLS analysis")
  cat("\n")
  cat("======================")
  cat("\n")

    cat("Global significance based on permuting the Escoufier RV coefficient")
    cat("\n")
    cat(paste0("Observed RV=", round(object$global_significance_RV["observed_RV"], 2)))
    cat("\n")
    cat(paste0("p value= ", object$global_significance_RV["RV_p_value"]))
    cat("\n")
    cat("Note that the observed RV is not rarefied")
    cat("\n")
    cat("======================")
    cat("\n")
    cat("Significance of each singular axis")
    cat("\n")
    SA_sign=object$singular_axis_significance
    SA_sign$correlation_PLS_scores=round(SA_sign$correlation_PLS_scores, 2)
    SA_sign$percentage_squared_covariance=round(SA_sign$percentage_squared_covariance, 2)
    if (all(SA_sign$singular_values>0.001)) {
      SA_sign$singular_values=round(SA_sign$singular_values, 3)
    }
    print(SA_sign)
    cat("======================")
  } else {
    print(object)
  }
}


# Function to perform the PLS fit
pls_base = function (x,y) {

  xor=x
  yor=y
  # original data

  x=scale(x, center = TRUE, scale = FALSE)
  y=scale(y, center = TRUE, scale = FALSE)
  # centered original data

  S12 = cov(x, y)
  # Obtain only the between-block covariance matrix

  SVD=corpcor::fast.svd(S12)
  # Use the fast.svd algorithm on the between-block covariance matrix

  XScores <- as.matrix(x) %*% SVD$u
  YScores <- as.matrix(y) %*% SVD$v
  # Compute scores for x and y


  CorrXScoresYScores=unlist(lapply(seq_len(length(SVD$d)), function(x) cor(XScores[,x], YScores[,x])))

  PLS=list(XScores = XScores, YScores = YScores,
           U=SVD$u, V=SVD$v, D=SVD$d,
           percentage_squared_covariance=(SVD$d^2/sum(SVD$d^2))*100,
           CorrXScoresYScores=CorrXScoresYScores,
           OriginalData=list(x=cbind(xor), y=cbind(yor)),
           x_center=attributes(x)$"scaled:center",
           y_center=attributes(y)$"scaled:center")
  class(PLS)=c("pls_fit", "list")
return(PLS)
}


# Function to run permutation tests
pls_perm = function (x, y, perm=999) {
  ObsPLS=pls_base(x,y)
  # Observed PLS analysis

  nobs=nrow(x)
  nsingval=length(ObsPLS$D)
  nvarX=ncol(x)
  nvarY=ncol(y)
  # Number of observations and singular values

    Yperm=lapply(seq_len(perm), function(k) y[sample(seq_len(nobs),nobs),])
  # Generate permuted samples for one of the two blocks
  # (the second)

  PLSpermD=matrix(rep(NA,nsingval*perm),nsingval,perm)
  PLSpermCorr=matrix(rep(NA,nsingval*perm),nsingval,perm)

  for (i in seq_len(perm)) {
    Temp=pls_base(x,Yperm[[i]])
    PLSpermD[,i]=Temp$D
    PLSpermCorr[,i]=Temp$CorrXScoresYScores
  }
  # Compute PLS on each of the permuted datasets and retain the singular values
  # Combine all permuted sets of singular values in columns
  # (so that each column is a different permutation,
  # each row a different singular value)

  D_p.values=rep(NA,nsingval)
  for (i in seq_len(nsingval)) {
    D_p.values[i]=(length(which(PLSpermD[i,]>=ObsPLS$D[i]))+1)/(perm+1)
  }
  # Compute p values as the proportion of permuted singular values larger
  # or equal to the one observed
  CorrXScoresYScores_p.values=rep(NA,nsingval)
  for (i in seq_len(nsingval)) {
    CorrXScoresYScores_p.values[i]=(length(which(PLSpermCorr[i,]>=ObsPLS$CorrXScoresYScores[i]))+1)/(perm+1)
  }
  # Same for the correlation between X scores and Y scores

  ObsRV=EscoufierRV(x, y)
  permRV=unlist(lapply(Yperm, function(Y) EscoufierRV(cbind(x),cbind(Y))))
  RV_p_value=(length(which(permRV>=ObsRV))+1)/(perm+1)
  # Overall significance using Escoufier RV

  Results=list(XScores = ObsPLS$XScores,
               YScores = ObsPLS$YScores,
               U=ObsPLS$U, V=ObsPLS$V, D=ObsPLS$D,
               CorrXScoresYScores=ObsPLS$CorrXScoresYScores,
               global_significance_RV=c(observed_RV=ObsRV, RV_p_value=RV_p_value),
               singular_axis_significance=data.frame(correlation_PLS_scores=ObsPLS$CorrXScoresYScores,
                                       pvalue_correlation_PLS_scores=CorrXScoresYScores_p.values,
                                       singular_values=ObsPLS$D,
                                       pvalue_singular_values=D_p.values,
                                       percentage_squared_covariance=(ObsPLS$D^2/sum(ObsPLS$D^2))*100
                                       ),
               OriginalData=list(x=ObsPLS$OriginalData$x, y=ObsPLS$OriginalData$y),
               x_center=ObsPLS$x_center,
               y_center=ObsPLS$y_center
               )
  class(Results)=c("pls_fit", "list")
  return(Results)
}


#' Major axis predictions for partial least squares (PLS) analysis
#'
#' Project data on the major axis of PLS scores and
#' obtain associated predictions
#'
#'
#' This function acts on a pls_fit object obtained from the function pls.
#' More in detail, the function:
#'  \itemize{
#'   \item{Projects the original data onto the major axis for each pair of PLS axes (obtaining for each observation of the original data a score along this axis)}
#'   \item{For each observation (specimen) of the original data, obtains the shape predicted by its score along the major axis}
#'   \item{(optionally) If new data is provided, these data is first projected in the original PLS space and then the two operations above are performed on the new data}
#' }
#' A more in-depth explanation with a figure which allows for a more intuitive understanding
#' is provided in Fruciano et al under review (kindly contact Carmelo Fruciano to receive the figure and relevant text prior to publication)
#' The idea is to obtain individual-level estimates of the shape predicted by a PLS model.
#' This can be useful, for instance, to quantify to which extent the shape of a given individual from one group
#' resembles the shape that individual would have according to the model computed in another group.
#' This can be done by obtaining predictions with this function and then computing the distance between
#' the actual shape observed for each individual and its prediction obtained from this function.
#' This is, indeed, how this approach has been used in Fruciano et al (under review).
#'
#' The function returns a list with two or three main elements which are themselves lists.
#' The most useful elements for the final user are highlighted in boldface.
#'
#' \emph{original_major_axis_projection} is a list containing as many elements as specified in axes_to_use (default 1).
#' Each of this elements contains the details of the computation of the major axis
#' (as a PCA of PLS scores for a pair of axes), and in particular:
#' \itemize{
#'   \item{major_axis_rotation: }{eigenvector}
#'   \item{mean_pls_scores: }{the mean scores for that axis pair used in the computation}
#'   \item{pls_scale: }{the scaling factor used}
#'   \item{\strong{original_data_PLS_projection}: }{scores of the original data on the major axis}
#' }
#'
#' \emph{original_major_axis_predictions_reversed} contains the predictions of the PLS model for the original data,
#' back-transformed to the original space (i.e., if the original data was shape, this will be shape).
#' If axes_to_use > 1, these predictions will be based on the major axis computed for all pairs of axes considered.
#' This element has two sub-elements:
#' \itemize{
#' \item{\strong{Block1}: }{prediction for block 1}
#' \item{\strong{Block2}: }{prediction for block 2}
#' }
#'
#' \emph{new_data_results} is only returned when new data is provided and contains the results of the
#' analyses obtained using a previous PLS model on new data and, in particular:
#' \itemize{
#' \item{\strong{new_data_Xscores}: }{PLS scores of the new data using the old model for the first block}
#' \item{\strong{new_data_Yscores}: }{PLS scores of the new data using the old model for the second block}
#' \item{\strong{new_data_major_axis_proj}: }{Scores of the new data on the major axis computed
#'  using the PLS model provided in pls_object. If axes_to_use > 1, each column correspond to a separate major axis}
#'  \item{\strong{new_data_Block1_proj_prediction_revert}: }{Predictions for the Block1 of the new data
#'  obtained by first computing the major axis projections for the new data (as found in element new_data_major_axis_proj)
#'   and then back-transforming these projection to the original space (e.g., shape)}
#'  \item{\strong{new_data_Block2_proj_prediction_revert}: }{Predictions for the Block2 of the new data
#'  obtained by first computing the major axis projections for the new data (as found in element new_data_major_axis_proj)
#'   and then back-transforming these projection to the original space (e.g., shape)}
#' }
#'
#'
#' @section Citation:
#' If you use this function, please cite Fruciano et al. under review
#' (contact Carmelo Fruciano for the most recent reference)
#'
#' @param pls_object object of class "pls_fit" obtained from the function pls
#' @param new_data_x,new_data_y (optional) matrices or data frames containing new data
#' @param axes_to_use number of pairs of PLS axes to use in the computation
#' (by default, this is performed only on the first axis)
#' @param scale_PLS logical indicating whether PLS scores for different blocks
#' should be scaled prior to computing the major axis
#'
#' @seealso \code{\link{pls}}
#'
#' @return The function outputs a list with the following elements (please, see the Details section for explanations on their sub-elements):
#'  \describe{
#'   \item{original_major_axis_projection}{For each PLS axis pair,
#'   results of the computation of major axis and projection of the original data on each axis}
#'   \item{original_major_axis_predictions_reversed}{Data obtained back-transforming the scores
#'   on the major axis into the original space (e.g., shape)}
#'   \item{new_data_results}{(only if new data has been provided) PLS scores for the new data,
#'   scores of the new data on the major axis, preditions for the new data back-transformed into the original space (e.g., shape)}
#' }
#'
#'
#'
#' @references Fruciano C, Colangelo P, Castiglia R, Franchini P. under review. Does divergence from normal patterns of integration increase as chromosomal fusions increase in number? A test on a house mouse hybrid zone.
#'
#' @section Notice:
#' \itemize{
#' \item{If new data is provided, this is first centered to the same average as in the original analysis, then it is translated back to the original scale}
#'}
#'
#' @examples
#'
#'
#'
#' ######################################
#' ### Example using the classical    ###
#' ### iris data set as a toy example ###
#' ######################################
#'
#' data(iris)
#' # Import the iris dataset
#' versicolor_data=iris[iris$Species=="versicolor",]
#' # Select only the specimens belonging to the species Iris versicolor
#' versicolor_sepal=versicolor_data[,grep("Sepal", colnames(versicolor_data))]
#' versicolor_petal=versicolor_data[,grep("Petal", colnames(versicolor_data))]
#' # Separate sepal and petal data for I. versicolor
#'
#'
#' PLS_sepal_petal_versicolor=pls(versicolor_sepal, versicolor_petal, perm=99)
#' summary(PLS_sepal_petal_versicolor)
#' # Compute the PLS for I. versicolor
#'
#'
#' plot(PLS_sepal_petal_versicolor$XScores[,1],
#'      PLS_sepal_petal_versicolor$YScores[,1],
#'      asp = 1,
#'      xlab = "PLS 1 Block 1 scores",
#'      ylab = "PLS 1 Block 2 scores")
#' # Plot the scores for the original data on the first pair of PLS axes (one axis per block)
#' # This is the data based on which we will compute the major axis direction
#' # Imagine fitting a line through those point, that is the major axis
#'
#' Pred_major_axis_versicolor=pls_major_axis(PLS_sepal_petal_versicolor, axes_to_use=1)
#' # Compute for I. versicolor the projections to the major axis
#' # using only the first pair of PLS axes (and scaling the scores prior to the computation)
#'
#' hist(Pred_major_axis_versicolor$original_major_axis_projection[[1]]$original_data_PLS_projection,
#'      main="Original data - projections on the major axis - based on the first pair of PLS axes",
#'      xlab="Major axis score")
#' # Plot distribution of PLS scores for each individual in the original data
#' # (I. versicolor)
#' # projected on the major axis for the first pair of PLS axis
#'
#' Pred_major_axis_versicolor$original_major_axis_predictions_reversed$Block1
#' Pred_major_axis_versicolor$original_major_axis_predictions_reversed$Block2
#' # Shape for each individual of the original data (I. versicolor)
#' # predicted by its position along the major axis
#'
#' # Now we will use the data from new species (I. setosa and I virginica)
#' # and obtain predictions from the PLS model obtained for I. versicolor
#'
#' # The easiest is to use the data for all three species
#' # as if they were both new data
#' # (using versicolor as new data is not going to affect the model)
#'
#'
#' all_sepal=iris[,grep("Sepal", colnames(iris))]
#' all_petal=iris[,grep("Petal", colnames(iris))]
#' # Separate sepals and petals (they are the two blocks)
#'
#' Pred_major_axis_versicolor_newdata=pls_major_axis(pls_object=PLS_sepal_petal_versicolor,
#'                                           new_data_x = all_sepal,
#'                                           new_data_y = all_sepal,
#'                                           axes_to_use=1)
#' # Perform the major axis computation using new data
#' # Notice that:
#' # - we are using the old PLS model (computed on versicolor only)
#' # - we are adding the new data in the same order as in the original model
#' #   (i.e., new_data_x is sepal data, new_data_y is petal data)
#'
#'
#' plot(Pred_major_axis_versicolor_newdata$new_data_results$new_data_Xscores[,1],
#'      Pred_major_axis_versicolor_newdata$new_data_results$new_data_Yscores[,1],
#'      col=iris$Species, asp=1,
#'      xlab = "Old PLS, Axis 1, Block 1",
#'      ylab = "Old PLS, Axis 1, Block 2")
#' # Plot the new data (both versicolor and setosa)
#' # in the space of the first pair of PLS axes computed only on versicolor
#' # The three species follow a quite similar trajectories
#' # but they have different average value on the major axis
#'
#' # To visualize this better, we can plot the scores along the major axis
#' # for the three species
#' boxplot(Pred_major_axis_versicolor_newdata$new_data_results$new_data_major_axis_proj[,1]~
#' iris$Species, xlab="Species", ylab="Score on the major axis")
#'
#' # We can also visualize the deviations from the major axis
#' # For instance by putting the predictions of the two blocks together
#' # Computing differences and then performing a PCA
#' predictions_all_data=cbind(
#'   Pred_major_axis_versicolor_newdata$new_data_results$new_data_Block1_proj_prediction_revert,
#'   Pred_major_axis_versicolor_newdata$new_data_results$new_data_Block2_proj_prediction_revert)
#' # Get the predictions for the two blocks (sepals and petals)
#' # and put them back together
#'
#' Euc_dist_from_predictions=unlist(lapply(seq(nrow(iris)), function(i)
#'   dist(rbind(iris[i,1:4],predictions_all_data[i,]))))
#' # for each flower, compute the Euclidean distance between
#' # the original values and what is predicted by the model
#'
#' boxplot(Euc_dist_from_predictions~iris$Species,
#'         xlab="Species", ylab="Euclidean distance from prediction")
#' # I. setosa is the one which deviates the most from the prediction
#'
#'
#'
#' @export
pls_major_axis=function(pls_object, new_data_x=NULL, new_data_y=NULL, axes_to_use=1, scale_PLS=TRUE) {
  if (!("pls_fit" %in% class(pls_object))) stop("pls_object should be an object of class 'pls_fit'")
  if (axes_to_use<1 | axes_to_use>ncol(pls_object$U)) stop("The axes_to_use should be comprised in the number of axes identified by the PLS analysis")
  # Initial checks

  PLS_scores_major_axis_find=lapply(seq(axes_to_use), function(i)
  {
    PLS_scores_temp=cbind(pls_object$XScores[,i],
                          pls_object$YScores[,i])
    if (scale_PLS==TRUE) {
      PCA_PLS_scores_temp=prcomp(PLS_scores_temp, scale. = TRUE)
    } else {
      PCA_PLS_scores_temp=prcomp(PLS_scores_temp, scale. = FALSE)
    }

    list(major_axis_rotation=cbind(PCA_PLS_scores_temp$rotation[,1]),
         mean_pls_scores=PCA_PLS_scores_temp$center,
         pls_scale=PCA_PLS_scores_temp$scale,
         original_data_PLS_projection=cbind(PCA_PLS_scores_temp$x[,1]))
  })
  names(PLS_scores_major_axis_find)=paste0("PLS_axis_pair_", seq(axes_to_use))

  # For each PLS axis in axes_to_use, perform a PCA of
  # PLS scores of one block and the other
  # to identify the major axis, as well as the projection
  # of the original PLS scores in each of these axes


  ### Now a double-reversion of the original data
  ### from projections on the major axis back to the original (shape) variables



  if (scale_PLS==TRUE) {
    revertPCAPLS=lapply(PLS_scores_major_axis_find, function(X) {
      Scores=cbind(X$original_data_PLS_projection)
      Eigenvectors=cbind(X$major_axis_rotation)
    ZEt = Scores %*% t(Eigenvectors)
    for (i in seq(ncol(ZEt))) {
      ZEt[,i]=ZEt[,i]*X$pls_scale[i]
    }
    rows = nrow(Scores)
    ZEt + rep.row(X$mean_pls_scores, rows)
    })
  } else {

  revertPCAPLS=lapply(PLS_scores_major_axis_find, function(X) {
    reversePCA(Scores = cbind(X$original_data_PLS_projection),
               Eigenvectors = cbind(X$major_axis_rotation),
               Mean = X$mean_pls_scores
    )
  })
  }
  # Get back scores in the PLS space for the major axis projections


  PLS_scores_major_axis_X=do.call("cbind", lapply(revertPCAPLS, function(X) X[,1]))
  PLS_scores_major_axis_Y=do.call("cbind", lapply(revertPCAPLS, function(X) X[,2]))
  # Separate X and Y (the two blocks) and put them together as column vectors

  Block1_revert=reversePCA(cbind(PLS_scores_major_axis_X),
                           Eigenvectors = cbind(pls_object$U[,seq(axes_to_use)]),
                           Mean = pls_object$x_center)
  Block2_revert=reversePCA(cbind(PLS_scores_major_axis_Y),
                           Eigenvectors = cbind(pls_object$V[,seq(axes_to_use)]),
                           Mean = pls_object$y_center)
  # Revert back to the original space (e.g., shape)

  Results=list(original_major_axis_projection=PLS_scores_major_axis_find,
               original_major_axis_predictions_reversed=list(Block1=Block1_revert,
                                                             Block2=Block2_revert)
  )

  # Now the analysis in case there is new data to project in the original
  # PLS

  if(!is.null(new_data_x) && !is.null(new_data_y)) {

    if ((ncol(new_data_x)!=ncol(pls_object$OriginalData$x)) |
        (ncol(new_data_y)!=ncol(pls_object$OriginalData$y)) ) {
      stop("The number of variables of the new data should be the same as the number of variables in the original data")
    }
    # Checks on consistency on the number of variables for new data

    new_data_x_center=as.matrix(sweep(cbind(new_data_x), 2, pls_object$x_center, "-"))
    new_data_y_center=as.matrix(sweep(cbind(new_data_y), 2, pls_object$y_center, "-"))
    # Translate the new data to the same origin as the original data


    new_data_Xscores=cbind(new_data_x_center %*% pls_object$U)
    new_data_Yscores=cbind(new_data_y_center %*% pls_object$V)
    # Compute PLS scores for the new data in the old PLS

    new_data_major_axis_proj=lapply(seq(axes_to_use), function(i)
    {
      PLS_scores_temp=cbind(new_data_Xscores[,i], new_data_Yscores[,i])
      new_data_mjp=PLS_scores_temp %*% PLS_scores_major_axis_find[[i]]$major_axis_rotation
    })
    names(new_data_major_axis_proj)=paste0("PLS_axis_pair_", seq(axes_to_use))
    # Compute projection of the PLS scores for the new data on the major axis of the old PLS

    new_data_revertPCAPLS=lapply(seq(axes_to_use), function(i) {
      reversePCA(Scores = new_data_major_axis_proj[[i]],
                 Eigenvectors = cbind(PLS_scores_major_axis_find[[i]]$major_axis_rotation),
                 Mean = PLS_scores_major_axis_find[[i]]$mean_pls_scores
      )
    })
    # Get back scores in the PLS space for the major axis projections of new data

    new_data_revertPCAPLS_major_axis_X=do.call("cbind", lapply(new_data_revertPCAPLS, function(X) X[,1]))
    new_data_revertPCAPLS_major_axis_Y=do.call("cbind", lapply(new_data_revertPCAPLS, function(X) X[,2]))
    # Separate X and Y (the two blocks) and put them together as column vectors

    new_data_Block1_revert=reversePCA(cbind(new_data_revertPCAPLS_major_axis_X),
                                      Eigenvectors = cbind(pls_object$U[,seq(axes_to_use)]),
                                      Mean = pls_object$x_center)
    new_data_Block2_revert=reversePCA(cbind(new_data_revertPCAPLS_major_axis_Y),
                                      Eigenvectors = cbind(pls_object$V[,seq(axes_to_use)]),
                                      Mean = pls_object$y_center)
    # Revert back to the original space (e.g., shape)

    Results$new_data_results=list(new_data_Xscores=new_data_Xscores,
                         new_data_Yscores=new_data_Yscores,
                         new_data_major_axis_proj=do.call("cbind", new_data_major_axis_proj),
                         new_data_Block1_proj_prediction_revert=new_data_Block1_revert,
                         new_data_Block2_proj_prediction_revert=new_data_Block2_revert
    )

  }


  return(Results)
}




