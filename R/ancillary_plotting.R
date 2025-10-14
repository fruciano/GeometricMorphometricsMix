### Ancillary plotting functions
### (not exported)

## Register internal helper column names to avoid R CMD check notes
utils::globalVariables(c("CI_point_color", "CI_errorbar_color", ".data"))

#' Create confidence interval plot
#'
#' Creates a simple ggplot showing points with error bars for
#' confidence intervals. Based on a data frame with group, average
#' values, and confidence interval limits.
#'
#' @param data A data frame containing the plotting data
#' @param x_var Character string specifying the column name for x-axis
#'   values (default "group")
#' @param y_var Character string specifying the column name for y-axis
#'   values (default "average")
#' @param ymin_var Character string specifying the column name for
#'   lower CI limit (default "CI_min")
#' @param ymax_var Character string specifying the column name for
#'   upper CI limit (default "CI_max")
#' @param x_lab Character string for x-axis label (default "Group")
#' @param y_lab Character string for y-axis label (default "Average")
#' @param point_color A single color or a vector of colors for point
#'   estimates. If length 1, the same color is used for all points.
#'   If length equals the number of unique x-axis levels OR the number
#'   of rows in `data`, colors are assigned per level (preferred) or
#'   per row respectively. (default "darkblue")
#' @param errorbar_color A single color or a vector of colors for
#'   error bars. Follows the same recycling / matching rules as
#'   `point_color`. (default "darkred")
#' @param ... Additional arguments passed to ggplot
#'
#' @return A ggplot object
#'
#' @noRd
CI_plot=function(data, x_var="group", y_var="average", 
                 ymin_var="CI_min", ymax_var="CI_max",
                 x_lab="Group", y_lab="Average", 
                 point_color="darkblue", errorbar_color="darkred", ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for CI_plot")
  }
  
  # Determine x levels for color assignment logic
  x_vals = data[[x_var]]
  x_levels = unique(x_vals)
  n_levels = length(x_levels)
  n_rows = nrow(data)

  # Helper to expand / validate a color vector
  expand_col_vector = function(col_vec, label) {
    if (length(col_vec) == 1) {
      return(rep(col_vec, n_rows))
    } else if (length(col_vec) == n_levels) {
      # Map by level
      named = setNames(col_vec, x_levels)
      return(unname(named[match(x_vals, names(named))]))
    } else if (length(col_vec) == n_rows) {
      return(col_vec)
    } else {
      stop(sprintf(
        "Length of %s (%d) must be 1, number of unique '%s' ",
        "levels (%d), or number of rows (%d)",
        label, length(col_vec), x_var, n_levels, n_rows
      ))
    }
  }

  # Prepare per-row colors (even if constant) for unified handling
  point_cols_row = expand_col_vector(point_color, "point_color")
  errorbar_cols_row = expand_col_vector(errorbar_color, "errorbar_color")

  data$CI_point_color = point_cols_row
  data$CI_errorbar_color = errorbar_cols_row

  # Build plot (no need for group=1)
  p = ggplot2::ggplot(
    data, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]]),
    ...
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = CI_point_color), alpha = 0.8, size = 3,
      show.legend = FALSE
    ) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data[[ymin_var]], ymax = .data[[ymax_var]], color = CI_errorbar_color), width = 0.1, show.legend = FALSE) +
    ggplot2::scale_color_identity() +
    ggplot2::theme_classic() +
    ggplot2::labs(x = x_lab, y = y_lab) +
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
