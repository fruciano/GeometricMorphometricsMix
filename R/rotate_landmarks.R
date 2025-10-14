#' User-defined rotation of a landmark configuration
#'
#' Rotates a single 2D or 3D landmark configuration about the
#' origin of the coordinate system.
#'
#' This function can be useful to change the orientation of data
#' for display purposes
#' It could also be used to empirically check the rotation-invariance
#' of an analysis.
#' Notice that this function works on a single configuration of
#' landmarks provided as a k x p matrix of k landmarks in p
#' dimensions (i.e., p is 2 for 2D data and 3 for 3D data)
#' As user-supplied rotation, the function expects a single number
#' in the case of 2D data (rotation on the plane),
#' a vector with three values (corresponding to rotation relative to
#' the three axes) in 3D.
#' If center=TRUE (default), first the configuration of landmarks is
#' centered, then the rotation is performed,
#' and finally the landmark coordinated are translated back to the
#' original position.
#' This accomplishes rotating the landmark configuration around its
#' center
#'
#' @param LMdata matrix k x p of k landmarks and p dimensions
#' @param rotation vector containing the rotation angle(s) (a single
#'   rotation for 2D data, three rotations for 3D data)
#' @param radians a logical on whether the angle(s) are provided in
#'   radians (default) or not
#' @param center a logical on whether the rotation should be on the
#'   centered configuration
#' @return The function outputs a matrix k x p of the original
#'   configuration of landmarks rotated according to the user-supplied
#'   rotation
#'
#'
#' @examples
#' library(ggplot2)
#'
#' Poly1=scale(t(matrix(c(4,1,3,1,1,2,2,3),nrow=2,ncol=4)),
			#' center=TRUE, scale=FALSE)
#' # Create a polygon centered at the origin
#' Poly2=rotate_landmarks(LMdata=Poly1, rotation=10, radians=FALSE)
#' # Create a new polygon which is the rotated version of the first
#' # with respect to the origin (rotation of 10 degrees)
#'
#' BothPolys4Plotting=as.data.frame(rbind(Poly1,Poly2))
#' BothPolys4Plotting[,3]=c(rep("Original",4),rep("Rotated",4))
#' BothPolys4Plotting[,4]=rep(1:4,2)
#' colnames(BothPolys4Plotting)=c("X","Y","Polygon","Landmark")
#' # Put them together in a way that is easy to plot in ggplot2
#' GraphLims=range(BothPolys4Plotting[,1:2])
#' # limits of the plot
#'
#' ggplot() +
#' geom_polygon(data=BothPolys4Plotting,
#'              mapping=aes(x=X, y=Y, group=Polygon, fill=Polygon),
#'              alpha=0.5) +
#' geom_point(data=BothPolys4Plotting, aes(x=X, y=Y, color=Polygon)) +
#' geom_text(data=BothPolys4Plotting, aes(x=X, y=Y, label=Landmark),
#'           hjust=1, vjust=1, size=4)+
#' coord_fixed(ratio=1, xlim=GraphLims, ylim=GraphLims)+
#' theme_classic()
#' # Plot the two polygons (landmarks are numbered for ease of
#' # visualization)
#'
#'
#' @import stats
#' @export
rotate_landmarks=function(LMdata, rotation, radians=TRUE, center=TRUE) {
	if (ncol(LMdata)==2 & length(rotation)!=1) {
	stop(paste("The dimensionality implied by landmark data and",
	           "rotation do not match"))
	}
	if (ncol(LMdata)==3 & length(rotation)!=3) {
	stop(paste("The dimensionality implied by landmark data and",
	           "rotation do not match"))
	}
	if (ncol(LMdata)>3 | ncol(LMdata)<2) {
	stop(paste("This function is valid for 2D or 3D data"))
	}
	# Error handling in case landmark data and/or rotation vector
	# do not match

	if (radians==FALSE) {
	rotation=vapply(rotation, deg2rad, FUN.VALUE = numeric(1))
	}

  if (center==TRUE) {
    linmod=lm(LMdata~1)
    LMdata=linmod$residuals
  }
  # Center the landmark configuration in case
  # a rotation around the centroid is desidered

	if (ncol(LMdata)==2) {
	cos_rot=cos(rotation)
	sin_rot=sin(rotation)
	rotmat=rbind(c(cos_rot, -sin_rot), c(sin_rot, cos_rot))
	# rotation matrix
	TransformedData=LMdata %*% rotmat
	}
	if (ncol(LMdata)==3) {
	Rx=rbind(c(1,0,0),
	         c(0, cos(rotation[1]), -sin(rotation[1])),
	         c(0, sin(rotation[1]), cos(rotation[1])))
	Ry=rbind(c(cos(rotation[2]), 0, sin(rotation[2])),
	         c(0,1,0),
	         c(-sin(rotation[2]), 0, cos(rotation[2])))
	Rz=rbind(c(cos(rotation[3]), -sin(rotation[3]), 0),
	         c(sin(rotation[3]), cos(rotation[3]), 0),
	         c(0,0,1))
	rotmat=Rz%*%Ry%*%Rx
	TransformedData=LMdata %*% rotmat
	}

  if (center==TRUE) {
    TransformedData=TransformedData+linmod$fitted.values
  }
  # If the rotation is around the centroid,
  # translate back to the original position

return(TransformedData)
}
