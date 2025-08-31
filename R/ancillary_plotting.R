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


#' Create grouped density plot
#'
#' Creates a ggplot showing density plots for different groups with separate colors.
#' This function provides a convenient way to visualize the distribution of a continuous
#' variable across different groups.
#'
#' @param data A data frame containing the plotting data
#' @param value_var Character string specifying the column name for the values to plot density for
#' @param group_var Character string specifying the column name for grouping variable
#' @param alpha Transparency level for the density plots (default 0.25)
#' @param x_lab Character string for x-axis label (default uses value_var name)
#' @param y_lab Character string for y-axis label (default "Density")
#' @param title Character string for plot title (default NULL)
#' @param ... Additional arguments passed to ggplot
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' # Create example data
#' data = data.frame(
#'   values = c(rnorm(100, 0, 1), rnorm(100, 2, 1)),
#'   group = rep(c("A", "B"), each = 100)
#' )
#' 
#' # Create density plot
#' density_plot_group(data, "values", "group")
#' 
#' # With custom labels and transparency
#' density_plot_group(data, "values", "group", 
#'                   alpha = 0.5, 
#'                   x_lab = "Measurement", 
#'                   title = "Distribution by Group")
#' }
#'
#' @noRd
density_plot_group = function(data, value_var, group_var, alpha = 0.25, 
                             x_lab = NULL, y_lab = "Density", title = NULL, ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for density_plot_group")
  }
  
  # Set default x-axis label if not provided
  if (is.null(x_lab)) {
    x_lab = value_var
  }
  
  # Create the plot
  p = ggplot2::ggplot(data, ggplot2::aes_string(x = value_var, fill = group_var), ...) +
    ggplot2::geom_density(alpha = alpha) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = x_lab, y = y_lab, fill = group_var, title = title) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  return(p)
}
