test_that("ProjectOrthogonal removes the requested vector subspace", {
  Data = rbind(c(1, 2, 3), c(4, -1, 2), c(-2, 5, 7))
  direction = c(1, 2, 0)

  projected = ProjectOrthogonal(Data, direction)

  expect_equal(drop(projected %*% direction), rep(0, nrow(Data)), tolerance = 1e-12)
  expect_equal(ProjectOrthogonal(projected, direction), projected, tolerance = 1e-12)

  two_directions = cbind(c(1, 0, 0), c(0, 1, 0))
  projected_to_z = ProjectOrthogonal(Data, two_directions)
  expect_equal(projected_to_z, cbind(rep(0, nrow(Data)), rep(0, nrow(Data)), Data[, 3]),
               tolerance = 1e-12)
})

test_that("ProjectOrthogonal rejects a vector with the wrong dimension", {
  expect_error(ProjectOrthogonal(matrix(1:6, nrow = 2), c(1, 2)),
               "number of columns of Data")
})

test_that("reversePCA reconstructs full scores and adds the mean for partial scores", {
  Mean = c(10, -2, 4)
  Eigenvectors = diag(3)
  Scores = rbind(c(1, 2, 3), c(-4, 0.5, 2))

  reconstructed = reversePCA(Scores, Eigenvectors, Mean)
  expect_equal(reconstructed, sweep(Scores, 2, Mean, `+`))

  partial = reversePCA(Scores[, 1, drop = FALSE], Eigenvectors[, 1, drop = FALSE], Mean)
  expect_equal(partial, cbind(Mean[1] + Scores[, 1], rep(Mean[2], nrow(Scores)),
                              rep(Mean[3], nrow(Scores))))

  rotated_eigenvectors = rbind(c(0, -1, 0), c(1, 0, 0), c(0, 0, 1))
  rotated_scores = rbind(c(2, 3, 4), c(-1, 5, 0.5))
  rotated_mean = c(10, 20, -2)
  rotated_expected = rbind(c(7, 22, 2), c(5, 19, -1.5))

  expect_equal(reversePCA(rotated_scores, rotated_eigenvectors, rotated_mean),
               rotated_expected)
})

test_that("rotate_landmarks preserves distances and rotates about the centroid", {
  landmarks = rbind(c(2, 3), c(4, 3), c(2, 6))
  centroid = colMeans(landmarks)

  rotated = rotate_landmarks(landmarks, rotation = pi / 2, center = TRUE)

  expect_equal(colMeans(rotated), centroid, tolerance = 1e-12)
  expect_equal(as.vector(dist(rotated)), as.vector(dist(landmarks)), tolerance = 1e-12)

  rotated_about_origin = rotate_landmarks(matrix(c(1, 0), nrow = 1),
                                           rotation = pi / 2, center = FALSE)
  expect_equal(rotated_about_origin, matrix(c(0, -1), nrow = 1), tolerance = 1e-12)
})

test_that("rotate_landmarks supports 3D rotations and degree input", {
  landmarks_3d = rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
  rotated_3d = rotate_landmarks(landmarks_3d, rotation = c(0, 0, pi / 2), center = FALSE)
  expected_3d = rbind(c(0, -1, 0), c(1, 0, 0), c(0, 0, 1))

  expect_equal(rotated_3d, expected_3d, tolerance = 1e-12)
  expect_equal(rotate_landmarks(matrix(c(1, 0), nrow = 1), rotation = 90,
                                radians = FALSE, center = FALSE),
               expected_3d[1, 1:2, drop = FALSE], tolerance = 1e-12)
})

test_that("rotate_landmarks reports incompatible and unsupported dimensions", {
  expect_error(rotate_landmarks(matrix(1:6, ncol = 3), rotation = 0),
               "dimensionality implied")
  expect_error(rotate_landmarks(matrix(1:8, ncol = 4), rotation = 0), "valid for 2D or 3D")
})

test_that("rescale_by_landmark_distance sets each specimen's requested distance", {
  Data = array(0, dim = c(3, 2, 2))
  Data[, , 1] = rbind(c(0, 0), c(3, 4), c(1, 1))
  Data[, , 2] = rbind(c(-1, 2), c(1, 2), c(0, 5))
  target_lengths = c(10, 1)

  rescaled = rescale_by_landmark_distance(Data, lm1 = 1, lm2 = 2,
                                           lengths = target_lengths)

  expect_equal(dim(rescaled), dim(Data))
  for (specimen in seq_len(dim(Data)[3])) {
    measured = sqrt(sum((rescaled[1, , specimen] - rescaled[2, , specimen])^2))
    expect_equal(measured, target_lengths[specimen], tolerance = 1e-12)
    expect_equal(rescaled[, , specimen],
                 Data[, , specimen] * (target_lengths[specimen] /
                                         sqrt(sum((Data[1, , specimen] - Data[2, , specimen])^2))),
                 tolerance = 1e-12)
  }
})

test_that("LM_relativepos_check identifies changed coordinate order with specimen names", {
  Dataset = array(NA_real_, dim = c(3, 2, 3),
                  dimnames = list(NULL, NULL, c("reference", "translated", "switched")))
  Dataset[, , 1] = rbind(c(0, 0), c(1, 1), c(2, 2))
  Dataset[, , 2] = Dataset[, , 1] + 5
  Dataset[, , 3] = rbind(c(0, 0), c(2, 1), c(1, 2))

  result = LM_relativepos_check(Dataset, LM_to_check = 1:3)

  expect_identical(result$potentially_problematic_idx, 3L)
  expect_named(result$potentially_problematic_LMs, "switched")
  expect_equal(result$potentially_problematic_LMs[[1]][, 1:2],
               matrix(c(2, 3), nrow = 1), ignore_attr = TRUE)
  expect_equal(unname(result$potentially_problematic_LMs[[1]][1, "Coord1"]), 0)
  expect_equal(unname(result$potentially_problematic_LMs[[1]][1, "Coord2"]), 1)
})

test_that("LM_relativepos_check reports a clean dataset and validates input dimensions", {
  Dataset = array(NA_real_, dim = c(2, 2, 2))
  Dataset[, , 1] = rbind(c(0, 0), c(1, 2))
  Dataset[, , 2] = Dataset[, , 1] + 3

  expect_message({
    result = LM_relativepos_check(Dataset, LM_to_check = 1:2)
  }, "Cannot find potentially problematic specimens")
  expect_null(result)
  expect_error(LM_relativepos_check(c(1, 2, 3), LM_to_check = 1:2),
               "Dataset should be a k x p x n array")
  expect_error(LM_relativepos_check(array(1:4, dim = c(2, 1, 2)), LM_to_check = 1:2),
               "only 2D and 3D")
})
