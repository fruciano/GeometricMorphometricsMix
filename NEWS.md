# GeometricMorphometricsMix 0.6.0.0

* Cleaned up code (e.g., line length)
* Fixed internal use of `sapply()`
* Removed previously defunct functions `rarefied_convex_hull()` and `rarefied_disparty()`

# GeometricMorphometricsMix 0.5.1.0

* Modified `RVcomparison()` to accept grouped data with pairwise comparisons
* Improved input validation for RV coefficient functions

# GeometricMorphometricsMix 0.5.0.1

* Fixed scaling issue in `pls_major_axis()`

# GeometricMorphometricsMix 0.5

* Enhanced `RVrarefied()` to support multiple groups and customizable confidence intervals
* Added plotting functions for resampling estimates
* Changed deprecated disparity functions to defunct
* Created unit tests for `RVrarefied()`

# GeometricMorphometricsMix 0.4.3.0

* Switched `Kmultparallel()` to use `future_lapply()` for improved parallelization
* Parallelized bootstrapping and rarefaction analyses

# GeometricMorphometricsMix 0.4.2.0

* Created general resampling function returning indexes
* Added bootstrap capability at lower sample sizes to `disparity_resample()`

# GeometricMorphometricsMix 0.4.1.0

* Added range reporting option to `disparity_resample()` when CI=1
* Removed input validation from main disparity function

# GeometricMorphometricsMix 0.4

* Added new datasets for examples and vignettes

# GeometricMorphometricsMix 0.3

* Added support for multiple datasets and tree sets to `Kmultparallel()`
* Migrated `Kmultparallel()` to use future for parallelism
* Added S3 plot method for `Kmultparallel()` output with density plots
* Removed parallelsugar from suggested packages

# GeometricMorphometricsMix 0.2

* Updated package description and title to reflect heterogeneous methods
* Introduced new disparity resampling function with S3 methods
* Deprecated old disparity resampling functions
* Created vignette for disparity analyses
