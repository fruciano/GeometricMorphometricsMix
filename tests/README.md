# Unit tests for GeometricMorphometricsMix

The suite uses testthat edition 3. New tests follow the package's `=` assignment convention and use withr to restore random-number settings and other temporary state.

## What is tested

- Existing suites: disparity_resample and its print/plot methods, RVrarefied and its print/plot methods, and RVcomparison.
- Geometry and transformations: ProjectOrthogonal, reversePCA, rotate_landmarks, rescale_by_landmark_distance, and LM_relativepos_check.
- Statistical calculations: pls and its print/summary methods, pls_major_axis, repeated_measures_test, scaled_variance_of_eigenvalues, and EscoufierRV.
- Resampling and comparisons: adjRand_test, BTailTest, critical_angle, TestOfAngle, disparity_test, dist_mean_boot, parallel_analysis and its summary method, and defunct Kmultparallel compatibility errors.
- Internal ancillary calculations: known distance, variance, replication, and angle-conversion answers.

Tests use small mathematical examples, independently derived expectations, reconstruction and invariance checks, reproducible permutations/bootstrap runs, output contracts, and selected invalid inputs. They do not establish statistical power, biological validity, or runtime speed. Some exceptional inputs remain untested or have unresolved defects; coverage is not a guarantee of correctness.

## Running tests

From the package root, with dependencies installed:

```r
testthat::test_local(".")
testthat::test_local(".", filter = "geometry-transformations")
```

Declared optional dependencies should be available for a full run. Individual tests may skip when their required optional dependencies are unavailable. Existing long-running tests may skip on CRAN.

## Measuring coverage

Install covr as a development tool (it is not required by package users). From the package root:

```r
coverage = covr::package_coverage(".", type = "tests")
covr::percent_coverage(coverage)
covr::tally_coverage(coverage, by = "line")
```

Use the same dependency versions, platform, and test settings for before/after comparisons. Tests-only coverage excludes vignette and example execution. Do not add exclusions merely to increase the percentage. Review both uncovered code and the quality of assertions on covered code. The revision working area contains the detailed before/after results and reproducible tooling, outside this package.
