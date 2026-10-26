# Small fixtures and replicate counts keep these tests fast.

disparity_compare_fixture=function() {
  X=rbind(c(0, 0), c(3, 4), c(6, 8))
  list(Data=rbind(X, 2*X),
       group=factor(rep(c("high", "low"), each=3),
                    levels=c("high", "low")),
       statistics=c("multivariate_variance",
                    "mean_pairwise_euclidean_distance"))
}

test_that("disparity measures and contrasts agree with known geometry", {
  withr::local_seed(11)
  f=disparity_compare_fixture()
  # Three points separated by 5, 5 and 10 units: variance=25, mean distance=20/3.
  expected=list(difference=c(-75, -20/3),
                absolute_difference=c(75, 20/3),
                ratio=c(1/4, 1/2), log_ratio=log(c(1/4, 1/2)))
  for (comparison in names(expected)) {
    result=disparity_compare(
      as.data.frame(f$Data), f$group, statistic=f$statistics,
      comparison=comparison, n_resamples=9
    )
    expect_equal(result$results$Observed_grp1, c(25, 20/3), tolerance=1e-10)
    expect_equal(result$results$Observed_grp2, c(100, 40/3), tolerance=1e-10)
    expect_equal(result$results$contrast, expected[[comparison]], tolerance=1e-10)
  }
})

test_that("group levels and missing membership define the focal comparison", {
  withr::local_seed(12)
  f=disparity_compare_fixture()
  reversed=disparity_compare(
    f$Data, factor(f$group, levels=c("low", "high")),
    comparison="log_ratio", n_resamples=9
  )
  expect_identical(reversed$group_levels, c("low", "high"))
  expect_equal(reversed$results$contrast, log(4), tolerance=1e-10)
  ratio=disparity_compare(f$Data, factor(f$group, levels=c("low", "high")),
                          comparison="ratio", n_resamples=9)
  expect_equal(ratio$results$contrast, 4)

  extra=rbind(f$Data, c(1000, -1000))
  membership=factor(c(as.character(f$group), NA),
                    levels=c("high", "low", "unused"))
  for (method in c("permutation", "bootstrap")) {
    focal=withr::with_seed(13, disparity_compare(
      f$Data, as.character(f$group), method=method,
      comparison="difference", n_resamples=9
    ))
    excluded=withr::with_seed(13, disparity_compare(
      extra, membership, method=method, comparison="difference", n_resamples=9
    ))
    expect_identical(excluded$group_levels, c("high", "low"))
    expect_identical(excluded$group_sizes, c(high=3L, low=3L))
    expect_identical(excluded$results, focal$results)
    expect_identical(excluded$replicate_contrasts, focal$replicate_contrasts)
  }
})

test_that("permutation output is reproducible and uses the documented tails", {
  withr::local_seed(14)
  f=disparity_compare_fixture()
  for (comparison in c("difference", "absolute_difference", "ratio", "log_ratio")) {
    run=function() disparity_compare(
      f$Data, f$group, statistic=f$statistics,
      comparison=comparison, n_resamples=9
    )
    first=withr::with_seed(15, run())
    second=withr::with_seed(15, run())
    expect_identical(first$results, second$results)
    expect_identical(first$replicate_contrasts, second$replicate_contrasts)
    expect_identical(dim(first$replicate_contrasts), c(9L, 2L))
    expect_identical(colnames(first$replicate_contrasts), f$statistics)
    for (j in seq_along(f$statistics)) {
      draws=first$replicate_contrasts[,j]
      observed=first$results$contrast[j]
      upper=(1+sum(draws >= observed))/10
      lower=(1+sum(draws <= observed))/10
      expected=if (comparison == "absolute_difference") upper else min(1, 2*min(lower, upper))
      expect_equal(first$results$p_value[j], expected)
    }
  }
  tied=disparity_compare(matrix(0, nrow=6, ncol=2), f$group,
                         comparison="absolute_difference", n_resamples=9)
  expect_equal(tied$results$p_value, 1)
})

test_that("bootstrap intervals match independent draws and handle undefined logs", {
  f=disparity_compare_fixture()
  # Independent formulas; do not call the package's internal disparity helpers.
  measures=function(X) {
    centered=sweep(X, 2, colMeans(X))
    pairs=utils::combn(seq_len(nrow(X)), 2)
    distances=sqrt(colSums(t(X[pairs[1,],,drop=FALSE]-X[pairs[2,],,drop=FALSE])^2))
    c(sum(centered^2)/(nrow(X)-1), mean(distances))
  }
  reference=withr::with_seed(16, t(replicate(19, {
    A=f$Data[sample(seq_len(3), 3, replace=TRUE),,drop=FALSE]
    B=f$Data[3+sample(seq_len(3), 3, replace=TRUE),,drop=FALSE]
    measures(A)-measures(B)
  })))
  result=withr::with_seed(16, disparity_compare(
    f$Data, f$group, statistic=f$statistics, comparison="difference",
    method="bootstrap", n_resamples=19
  ))
  expect_equal(unname(result$replicate_contrasts), reference, tolerance=1e-10)
  for (j in seq_along(f$statistics)) {
    bounds=stats::quantile(reference[,j], c(0.025, 0.5, 0.975), names=FALSE)
    actual=unname(unlist(result$results[j,c("CI_lower", "CI_median", "CI_upper")]))
    expect_equal(actual, bounds, tolerance=1e-10)
  }
  expect_true(all(is.na(result$results$p_value)))
  expect_equal(result$results$n_valid, c(19, 19))

  small=matrix(c(0, 1, 0, 2), ncol=1)
  expect_warning({
    failed=withr::with_seed(17, disparity_compare(
      small, rep(c("high", "low"), each=2),
      method="bootstrap", comparison="log_ratio", n_resamples=19
    ))
  }, regexp="Undefined replicate contrasts")
  expect_true(is.finite(failed$results$contrast))
  expect_equal(failed$results$n_valid, sum(is.finite(failed$replicate_contrasts)))
  expect_lt(failed$results$n_valid, 19)
  expect_true(anyNA(failed$replicate_contrasts))
  expect_true(all(is.na(unlist(failed$results[c(
    "bootstrap_mean", "bootstrap_sd", "bootstrap_bias",
    "CI_lower", "CI_median", "CI_upper", "interval_contains_null"
  )]))))
})

test_that("invalid measurements, groups and resampling settings are rejected", {
  f=disparity_compare_fixture()
  expect_error(disparity_compare(replace(f$Data, 1, NA), f$group), "finite")
  expect_error(disparity_compare(f$Data, rep(c("a", "b", "c"), each=2)),
               "exactly two")
  expect_error(disparity_compare(f$Data, c("a", rep("b", 5))), "at least two")
  expect_error(disparity_compare(f$Data, f$group[-1]), "one entry per row")
  for (n in c(0, -1, 1.5, NA_real_)) {
    expect_error(disparity_compare(f$Data, f$group, n_resamples=n), "positive integer")
  }
  expect_error(disparity_compare(f$Data, f$group, method="bootstrap", n_resamples=1),
               "at least two replicates")
  expect_error(disparity_compare(matrix(0, nrow=6, ncol=2), f$group,
                                 comparison="log_ratio"), "Undefined observed contrast")
})

test_that("S3 printing is compact and leaves results and RNG unchanged", {
  withr::local_seed(18)
  f=disparity_compare_fixture()
  expect_true(is.function(getS3method("print", "disparity_compare")))
  for (method in c("permutation", "bootstrap")) {
    result=disparity_compare(f$Data, f$group, statistic=f$statistics,
                             comparison="difference", method=method, n_resamples=9)
    expect_s3_class(result, "disparity_compare")
    before=result
    seed_before=.Random.seed
    expect_warning({
      output=capture.output({
        visible=withVisible(print(result))
      })
    }, regexp=NA)
    expect_false(visible$visible)
    expect_identical(visible$value, before)
    expect_identical(result, before)
    expect_identical(.Random.seed, seed_before)
    expect_match(paste(output, collapse="\n"), "Disparity comparison", fixed=TRUE)
    expect_match(paste(output, collapse="\n"), "Valid/Requested", fixed=TRUE)
    expect_false(any(grepl("replicate_contrasts|\\$fit|\\$call", output)))
  }
  expect_warning({
    failed=withr::with_seed(17, disparity_compare(
      matrix(c(0, 1, 0, 2), ncol=1), rep(c("high", "low"), each=2),
      method="bootstrap", comparison="log_ratio", n_resamples=19
    ))
  }, regexp="Undefined replicate contrasts")
  expect_warning({
    output=capture.output(print(failed))
  }, regexp=NA)
  expect_true(any(grepl("NA", output, fixed=TRUE)))
  expect_false(any(grepl("withheld|Undefined", output)))
})
