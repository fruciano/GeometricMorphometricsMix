test_that("internal distance and variance helpers agree with elementary answers", {
  x = rbind(c(0, 0), c(3, 0), c(0, 4))
  expect_equal(GeometricMorphometricsMix:::muvar(x), sum(diag(stats::cov(x))))
  expect_equal(GeometricMorphometricsMix:::meanpairwiseEuclideanD(x), 4)
  expect_equal(as.numeric(GeometricMorphometricsMix:::pdistance(c(0, 0), c(3, 4))), 5)
  expect_equal(GeometricMorphometricsMix:::rad2deg(pi), 180)
  expect_equal(GeometricMorphometricsMix:::deg2rad(90), pi / 2)
  expect_equal(GeometricMorphometricsMix:::rep_row(c(2, 5), 3),
               matrix(c(2, 2, 2, 5, 5, 5), nrow = 3))
  expect_equal(GeometricMorphometricsMix:::rep_col(c(2, 5), 3),
               matrix(c(2, 5, 2, 5, 2, 5), nrow = 2))
})
