library(testthat)
library(ape)
library(GeometricMorphometricsMix)

# Ensure sequential future plan in tests
if (requireNamespace("future", quietly = TRUE)) future::plan(future::sequential)

test_that("single dataset and treeset returns expected columns", {
  set.seed(1)
  tree <- rtree(10)
  # wrap single phylo into a multiPhylo
  tree <- list(tree)
  class(tree) <- "multiPhylo"
  data_mat <- matrix(rnorm(10 * 3), nrow = 10, ncol = 3)
  rownames(data_mat) <- tree[[1]]$tip.label

  res <- Kmultparallel(data_mat, tree, verbose = FALSE)
  expect_s3_class(res, c("parallel_Kmult", "data.frame"))
  expect_true(all(c("Kmult", "treeset", "dataset", "tree_index") %in% colnames(res)))
  expect_true(nrow(res) == length(tree))
})


test_that("list of datasets and treesets works and dimensions match", {
  set.seed(2)
  tree1 <- rtree(8); tree2 <- rtree(8)
  tree1$tip.label <- paste0("sp", 1:8)
  tree2$tip.label <- paste0("sp", 1:8)
  multi1 <- list(tree1); class(multi1) <- "multiPhylo"
  multi2 <- list(tree2); class(multi2) <- "multiPhylo"
  treeset <- list(ts1 = multi1, ts2 = multi2)

  data1 <- matrix(rnorm(8*2), nrow = 8); rownames(data1) <- paste0("sp", 1:8)
  data2 <- matrix(rnorm(8*2), nrow = 8); rownames(data2) <- paste0("sp", 1:8)
  datasets <- list(d1 = data1, d2 = data2)

  res <- Kmultparallel(datasets, treeset, verbose = FALSE)
  # expect two datasets * two treeset entries (each treeset has 1 tree)
  expect_true(nrow(res) == 4)
  expect_true(all(res$dataset %in% names(datasets)))
  expect_true(all(res$treeset %in% names(treeset)))
})


test_that("error when data has fewer than 3 rows", {
  tree <- rtree(5)
  tree <- list(tree); class(tree) <- "multiPhylo"
  data_mat <- matrix(rnorm(2*3), nrow = 2)
  rownames(data_mat) <- tree[[1]]$tip.label[1:2]
  expect_error(Kmultparallel(data_mat, tree, verbose = FALSE), "fewer than 3 rows")
})


test_that("iter > 0 returns p value column with values in [0,1]", {
  set.seed(3)
  tree <- rtree(7)
  tree <- list(tree); class(tree) <- "multiPhylo"
  data_mat <- matrix(rnorm(7*2), nrow = 7); rownames(data_mat) <- tree[[1]]$tip.label
  res <- Kmultparallel(data_mat, tree, iter = 10, verbose = FALSE)
  expect_true("p value" %in% colnames(res))
  expect_true(all(res$`p value` >= 0 & res$`p value` <= 1))
})
