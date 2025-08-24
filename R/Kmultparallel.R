#' Parallel implementation of Adams' Kmult with additional support for multiple datasets and tree sets
#'
#' Parallel implementation of Kmult, a measure of phylogenetic signal which
#' is a multivariate equivalent of Blomberg's K. This version supports multiple
#' datasets and tree sets, computing Kmult for all combinations.
#'
#' This is an updated and improved version of the function included in Fruciano et al. 2017.
#' It performs the computation of Adams' Kmult (Adams 2014) in parallel
#' with the aim of facilitating computation on a distribution of trees rather than a single tree.
#' This version uses the future framework for parallelization, making it more flexible
#' and compatible across different operating systems.
#' If one wanted to perform a computation of Kmult on a single tree, he/she would be
#' advised to use the version implemented in the package geomorph, which receives regular updates.
#'
#' @section Parallelization:
#' This function uses the future framework for parallelization. To enable parallel processing,
#' you need to set up a future plan before calling this function. For example:
#' \itemize{
#'   \item For sequential processing: \code{future::plan(future::sequential)}
#'   \item For multicore processing (Unix/Mac): \code{future::plan(future::multicore, workers = 4)}
#'   \item For multisession processing (Windows/Unix/Mac): \code{future::plan(future::multisession, workers = 4)}
#'   \item For cluster processing: \code{future::plan(future::cluster, workers = c("host1", "host2"))}
#' }
#' On Windows machines, use \code{future::plan(future::multisession, workers = 4)} for parallel processing.
#' The number of workers should typically not exceed the number of CPU cores available.
#'
#' @section Citation:
#' If you use this function please kindly cite both
#' Fruciano et al. 2017 (because you're using this parallelized function) and
#' Adams 2014 (because the function computes Adams' Kmult)
#'
#' @param data Either a data.frame/matrix with continuous (multivariate) phenotypes,
#' or a list where each element is a data.frame/matrix representing a separate dataset.
#' Row names should match species names in the phylogenetic trees.
#' @param trees Either a multiPhylo object containing a collection of trees (single tree set),
#' or a list where each element is a multiPhylo object representing a separate tree set.
#' @param burninpercent percentage of trees in each tree set to discard as burn-in
#' (by default no tree is discarded)
#' @param iter number of permutations to be used in the permutation test
#' (this should normally be left at the default value of 0 as permutations slow down
#' computation and are of doubtful utility when analyzing tree distributions)
#'
#' @return The function outputs a data.frame with classes "parallel_Kmult" and "data.frame" containing columns:
#'  \describe{
#'   \item{Kmult}{Value of Kmult for each tree-dataset combination}
#'   \item{p value}{p value for the significance of the test (only if iter > 0)}
#'   \item{treeset}{Identifier for the tree set (name from list or number)}
#'   \item{dataset}{Identifier for the dataset (name from list or number)}
#'   \item{tree_index}{Index of the tree within its tree set}
#' }
#'
#' @section S3 Methods:
#' The returned object has specialized S3 methods:
#' \itemize{
#'   \item \code{\link{print.parallel_Kmult}}: Provides a summary of Kmult ranges for each dataset-treeset combination
#'   \item \code{\link{plot.parallel_Kmult}}: Creates density plots of Kmult values grouped by dataset-treeset combinations
#' }
#'
#' @examples
#' \dontrun{
#' # Set up future for parallel processing on multiple cores (including Windows)
#' library(future)
#' plan(multisession, workers = 6)
#' # Use multisession backend which works on all platforms including Windows
#' 
#' # Load required packages for data simulation
#' library(phytools)
#' library(MASS)
#' library(mvMORPH)
#' 
#' # Generate 20 random phylogenetic trees with 100 tips each
#' all_trees = replicate(20, pbtree(n = 100), simplify = FALSE)
#' class(all_trees) = "multiPhylo"
#' # Create a collection of 20 random trees
#' 
#' # Split trees into 2 tree sets
#' treeset1 = all_trees[1:5]
#' treeset2 = all_trees[6:20]
#' class(treeset1) = class(treeset2) = "multiPhylo"
#' # Split the 20 trees into 2 separate tree sets
#' 
#' # Get tip names from the first tree for consistent naming
#' tip_names = all_trees[[1]]$tip.label[1:40]
#' # Use first 40 tip names for consistent data generation
#' 
#' # Generate 2 datasets using multivariate normal distribution
#' dataset1 = mvrnorm(n = 40, mu = rep(0, 5), Sigma = diag(5))
#' dataset2 = mvrnorm(n = 40, mu = rep(1, 5), Sigma = diag(5) * 2)
#' rownames(dataset1) = rownames(dataset2) = tip_names
#' # Create two datasets with different means and variances
#' # Notice how these datasets are random and should not display any phylogenetic signal
#' 
#' # Generate 5 datasets using Brownian motion evolution on the first treeset
#' bm_datasets = lapply(1:5, function(i) {
#'   tree_temp = treeset1[[i]]
#'   # Get only the first 40 tips to match our data size
#'   tips_to_keep = tree_temp$tip.label[1:40]
#'   tree_pruned = drop.tip(tree_temp, setdiff(tree_temp$tip.label, tips_to_keep))
#'   
#'   # Simulate data under Brownian motion
#'   sim_data = mvSIM(tree = tree_pruned, nsim = 1, model = "BM1", 
#'                    param = list(sigma = diag(5), theta = rep(0, 5)))
#'   return(sim_data)
#' })
#' names(bm_datasets) = paste0("BM_dataset_", 1:5)
#' # Generate 5 datasets evolving under Brownian motion
#' # Notice that these datasets should display strong phylogenetic signal
#' # when each of them is combined with the individual tree used in the simulation of trait evolution
#' 
#' # Example 1: Single dataset and single treeset analysis
#' result_single = Kmultparallel(bm_datasets[[1]], treeset1)
#' # Analyze first BM dataset with first treeset
#' 
#' # Use S3 methods to examine results
#' print(result_single)
#' # Display summary of Kmult values
#' # Notice how the range is very broad because we have high phylogenetic signal for those cases
#' # in which a given dataset has been simulated under Brownian motion with a given tree, but low phylogenetic signal
#' # when we use a different tree.
#' 
#' plot(result_single)
#' # Create density plot of Kmult distribution
#' # Notice the bimodal distribution with low phylogenetic signal corresponding to a
#' # mismatch between the tree used and the true evolutionary history of the traits,
#' # and the high phylogenetic signal when the correct tree is used.
#' 
#' # Example 2: Multiple datasets and multiple treesets analysis
#' # Combine all datasets into a list
#' all_datasets = c(list(mvrnorm1 = dataset1, mvrnorm2 = dataset2), bm_datasets)
#' # Combine random and BM datasets
#' 
#' # Combine all treesets into a list
#' all_treesets = list(treeset1 = treeset1, treeset2 = treeset2, 
#'                     treeset3 = treeset3, treeset4 = treeset4)
#' # Create list of all tree sets
#' 
#' # Run comprehensive analysis on all combinations
#' result_multiple = Kmultparallel(all_datasets, all_treesets)
#' # Analyze all dataset-treeset combinations
#' 
#' # Examine results using S3 methods
#' print(result_multiple)
#' # Display summary showing ranges for each combination
#' 
#' plot(result_multiple)
#' # Create grouped density plots by combination
#' 
#' # Custom plotting with different transparency
#' plot(result_multiple, alpha = 0.5, title = "Kmult Distribution Across All Combinations")
#' # Customize the plot appearance
#' }
#'
#' @references Adams DC. 2014. A Generalized K Statistic for Estimating Phylogenetic Signal from Shape and Other High-Dimensional Multivariate Data. Systematic Biology 63:685-697.
#' @references Fruciano C, Celik MA, Butler K, Dooley T, Weisbecker V, Phillips MJ. 2017. Sharing is caring? Measurement error and the issues arising from combining 3D morphometric datasets. Ecology and Evolution 7:7034-7046.
#'
#'
#' @import stats
#' @importFrom ape drop.tip vcv.phylo
#' @importFrom future.apply future_lapply
#' @export
Kmultparallel = function(data, trees, burninpercent = 0, iter = 0) {
    
    # Standardize input formats
    # Convert single datasets/tree sets to lists for uniform processing
    if (is.data.frame(data) || is.matrix(data)) {
        data = list(dataset1 = data)
        data_names = "1"
    } else if (is.list(data)) {
        data_names = if (is.null(names(data))) as.character(seq_along(data)) else names(data)
    } else {
        stop("data must be a data.frame, matrix, or list of data.frames/matrices")
    }
    
    if (inherits(trees, "multiPhylo")) {
        trees = list(treeset1 = trees)
        tree_names = "1"
    } else if (is.list(trees)) {
        tree_names = if (is.null(names(trees))) as.character(seq_along(trees)) else names(trees)
        # Check that all elements are multiPhylo
        if (!all(sapply(trees, function(x) inherits(x, "multiPhylo")))) {
            stop("All elements in trees list must be multiPhylo objects")
        }
    } else {
        stop("trees must be a multiPhylo object or list of multiPhylo objects")
    }
    
    # Apply burn-in to each tree set
    trees_processed = lapply(trees, function(treeset) {
        burnin_start = round((length(treeset)/100) * burninpercent)
        if (burnin_start >= length(treeset)) burnin_start = 0
        treeset[(burnin_start + 1):length(treeset)]
    })
    
    # Create all combinations of datasets and tree sets
    combinations = expand.grid(
        dataset_idx = seq_along(data),
        treeset_idx = seq_along(trees_processed),
        stringsAsFactors = FALSE
    )
    
    # Function to process a single combination
    process_combination = function(combo_idx) {
        combo = combinations[combo_idx, ]
        current_data = data[[combo$dataset_idx]]
        current_trees = trees_processed[[combo$treeset_idx]]
        dataset_name = data_names[combo$dataset_idx]
        treeset_name = tree_names[combo$treeset_idx]
        
        # Check if all trees in the treeset have the same tips
        # Extract tip labels from all trees
        all_tip_labels = lapply(current_trees, function(tree) tree$tip.label)
        
        # Check if all tip sets are identical
        first_tips = all_tip_labels[[1]]
        all_same_tips = all(sapply(all_tip_labels, function(tips) {
            length(tips) == length(first_tips) && all(sort(tips) == sort(first_tips))
        }))
        
        if (all_same_tips) {
            # All trees have the same tips - efficient processing
            # Find tips to drop based on the first tree (same for all)
            droplist = setdiff(current_trees[[1]]$tip.label, row.names(current_data))
            
            # Function to prune a single tree (same droplist for all)
            prune_tree = function(tree) {
                if (length(droplist) > 0) {
                    ape::drop.tip(tree, droplist)
                } else {
                    tree
                }
            }
            
            # Prune all trees in the set
            pruned_trees = future.apply::future_lapply(current_trees, prune_tree, future.seed = TRUE)
            class(pruned_trees) = "multiPhylo"
            
            # Reorder data to match tree tip order for each tree
            data_reordered = lapply(pruned_trees, function(tree) {
                current_data[tree$tip.label, , drop = FALSE]
            })
        } else {
            # Trees have different tips - individual processing
            # Process each tree individually with its own droplist
            tree_processing_results = future.apply::future_lapply(current_trees, function(tree) {
                # Find tips to drop for this specific tree
                droplist_tree = setdiff(tree$tip.label, row.names(current_data))
                
                # Prune this tree
                if (length(droplist_tree) > 0) {
                    pruned_tree = ape::drop.tip(tree, droplist_tree)
                } else {
                    pruned_tree = tree
                }
                
                # Reorder data for this tree
                data_reordered_tree = current_data[pruned_tree$tip.label, , drop = FALSE]
                
                list(pruned_tree = pruned_tree, data_reordered = data_reordered_tree)
            }, future.seed = TRUE)
            
            # Extract pruned trees and reordered data
            pruned_trees = lapply(tree_processing_results, function(x) x$pruned_tree)
            class(pruned_trees) = "multiPhylo"
            data_reordered = lapply(tree_processing_results, function(x) x$data_reordered)
        }
        
        # Compute Kmult for each tree
        tree_indices = seq_along(pruned_trees)
        kmult_results = future.apply::future_lapply(tree_indices, function(i) {
            result = Test_Kmult(data_reordered[[i]], pruned_trees[[i]], iter = iter)
            list(
                kmult = result$phy.signal,
                pvalue = result$pvalue,
                treeset = treeset_name,
                dataset = dataset_name,
                tree_index = i
            )
        }, future.seed = TRUE)
        
        return(kmult_results)
    }
    
    # Process all combinations
    all_results = future.apply::future_lapply(seq_len(nrow(combinations)), process_combination, future.seed = TRUE)
    
    # Flatten results and convert to data.frame
    flattened_results = do.call(c, all_results)
    
    result_df = data.frame(
        Kmult = sapply(flattened_results, function(x) x$kmult),
        treeset = sapply(flattened_results, function(x) x$treeset),
        dataset = sapply(flattened_results, function(x) x$dataset),
        tree_index = sapply(flattened_results, function(x) x$tree_index),
        stringsAsFactors = FALSE
    )
    
    # Add p-value column only if iter > 0
    if (iter > 0) {
        result_df$"p value" = sapply(flattened_results, function(x) x$pvalue)
    }
    
    # Add classes to the result
    class(result_df) = c("parallel_Kmult", "data.frame")
    
    return(result_df)
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
    # hist(K.val, 30, freq = TRUE, col = "gray", xlab = "Phylogenetic Signal")
    # arrows(K.obs, 50, K.obs, 5, length = 0.1, lwd = 2)
    return(list(phy.signal = K.obs, pvalue = P.val))
}


#' Print method for parallel_Kmult objects
#'
#' Provides a summary of Kmult analysis results showing the range of Kmult values
#' for each combination of dataset and treeset.
#'
#' @param x An object of class 'parallel_Kmult' produced by \code{\link{Kmultparallel}}
#' @param ... Additional arguments (currently not used)
#'
#' @return Invisibly returns the input object
#'
#' @examples
#' \dontrun{
#' # Assuming you have data and trees
#' result = Kmultparallel(data, trees)
#' print(result)  # or simply: result
#' }
#'
#' @export
print.parallel_Kmult = function(x, ...) {
    cat("Parallel Kmult Analysis Results\n")
    cat("===============================\n\n")
    
    # Create combination identifier
    x$combination = paste0("Dataset:", x$dataset, " - Treeset:", x$treeset)
    
    # Get unique combinations
    unique_combinations = unique(x$combination)
    
    cat("Summary by Dataset-Treeset combination:\n\n")
    
    for (combo in unique_combinations) {
        subset_data = x[x$combination == combo, ]
        n_trees = nrow(subset_data)
        kmult_range = range(subset_data$Kmult)
        kmult_mean = mean(subset_data$Kmult)
        kmult_sd = sd(subset_data$Kmult)
        
        cat(sprintf("%s\n", combo))
        cat(sprintf("  Trees analyzed: %d\n", n_trees))
        cat(sprintf("  Kmult range: %.4f - %.4f\n", kmult_range[1], kmult_range[2]))
        cat(sprintf("  Kmult mean (SD): %.4f (%.4f)\n\n", kmult_mean, kmult_sd))
    }
    
    # Check if p-values are present
    if ("p value" %in% colnames(x)) {
        cat("Note: Permutation p-values are included in the results\n")
    }
    
    invisible(x)
}


#' Plot method for parallel_Kmult objects
#'
#' Creates density plots of Kmult values, with separate densities for each
#' combination of dataset and treeset (if multiple combinations are present).
#'
#' @param x An object of class 'parallel_Kmult' produced by \code{\link{Kmultparallel}}
#' @param alpha Transparency level for density plots (default 0.25)
#' @param title Character string for plot title (default NULL for automatic title)
#' @param x_lab Character string for x-axis label (default "Kmult")
#' @param ... Additional arguments passed to the plotting function
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' # Assuming you have data and trees
#' result = Kmultparallel(data, trees)
#' plot(result)
#' 
#' # With custom settings
#' plot(result, alpha = 0.5, title = "Kmult Distribution")
#' }
#'
#' @export
plot.parallel_Kmult = function(x, alpha = 0.25, title = NULL, x_lab = "Kmult", ...) {
    
    # Create combination identifier for grouping
    x$combination = paste0("Dataset:", x$dataset, " - Treeset:", x$treeset)
    
    # Check if there are multiple combinations
    unique_combinations = unique(x$combination)
    n_combinations = length(unique_combinations)
    
    # Set default title if not provided
    if (is.null(title)) {
        if (n_combinations == 1) {
            title = "Kmult Distribution"
        } else {
            title = "Kmult Distribution by Dataset-Treeset Combination"
        }
    }
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plotting parallel_Kmult objects")
    }
    
    if (n_combinations == 1) {
        # Single combination - simple density plot
        p = ggplot2::ggplot(x, ggplot2::aes(x = Kmult)) +
            ggplot2::geom_density(alpha = alpha, fill = "steelblue") +
            ggplot2::theme_classic() +
            ggplot2::labs(x = x_lab, y = "Density", title = title) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    } else {
        # Multiple combinations - use grouped density plot
        p = ggplot2::ggplot(x, ggplot2::aes(x = Kmult, fill = combination)) +
            ggplot2::geom_density(alpha = alpha) +
            ggplot2::theme_classic() +
            ggplot2::labs(x = x_lab, y = "Density", fill = "Combination", title = title) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }
    
    return(p)
}
