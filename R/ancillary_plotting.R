### Ancillary plotting functions
### (not exported)

#' Create confidence interval plot
#'
#' Creates a simple ggplot showing points with error bars for confidence intervals.
#' Based on a data frame with group, observed values, and confidence interval limits.
#'
#' @param data A data frame containing the plotting data
#' @param x_var Character string specifying the column name for x-axis values (default "group")
#' @param y_var Character string specifying the column name for y-axis values (default "observed")
#' @param ymin_var Character string specifying the column name for lower CI limit (default "CI_min")
#' @param ymax_var Character string specifying the column name for upper CI limit (default "CI_max")
#' @param x_lab Character string for x-axis label (default "Group")
#' @param y_lab Character string for y-axis label (default "Observed")
#' @param ... Additional arguments passed to ggplot
#'
#' @return A ggplot object
#'
#' @noRd
CI_plot=function(data, x_var="group", y_var="observed", 
                 ymin_var="CI_min", ymax_var="CI_max",
                 x_lab="Group", y_lab="Observed", ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for CI_plot")
  }
  
  # Create the plot
  p=ggplot2::ggplot(data, ggplot2::aes_string(x=x_var, y=y_var, group=1), ...) +
    ggplot2::geom_point(alpha=0.8, size=3, colour="darkblue") +
    ggplot2::theme_classic() +
    ggplot2::geom_errorbar(width=.1, 
                          ggplot2::aes_string(ymin=ymin_var, ymax=ymax_var), 
                          colour="darkred") +
    ggplot2::labs(x=x_lab, y=y_lab) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  return(p)
}
